#!/usr/bin/python3
"""
Download a KEGG map PNG file ("2x resolution"), the associated KGML KO file, and the associated
KGML RN file for a global map. Rescale the objects in the KGML files to fit a 2x map.

Author: Samuel Miller
Contact: samuelmiller10@gmail.com
"""

import os
import re
import argparse
import urllib.request

from typing import Tuple

kegg_rest_api_get = 'http://rest.kegg.jp/get'

map_id_pattern = re.compile(r'map\d{5}')
# Overview maps have a certain range of IDs.
overview_map_id_pattern = re.compile(r'map\d{1}1[12]\d{2}')

def get_args() -> argparse.Namespace:
    """
    Get command line arguments.
    
    Returns
    =======
    argparse.Namespace
        Command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'map_id',
        help="ID of the map to set up, formatted as 'mapXXXXX', where 'X' is a digit."
    )
    parser.add_argument(
        '--map_filepath',
        default=None,
        help=(
            "Download filepath for the 2x resolution map PNG file. By default, without this "
            "option, the filepath is './mapXXXXX.png', where 'mapXXXXX' is the specified map ID."
        )
    )
    parser.add_argument(
        '--ko_kgml_filepath',
        default=None,
        help=(
            "Download filepath for the map KGML KO file. The file is modified in place to fit to a "
            "2x resolution map. By default, without this option, the filepath is './koXXXXX.xml', "
            "where 'koXXXXX' is derived from the specified map ID."
        )
    )
    parser.add_argument(
        '--rn_kgml_filepath',
        default=None,
        help=(
            "Download filepath for the map KGML RN file, used in drawing overview maps. The file "
            "is modified in place to fit to a 2x resolution map. By default, without this option, "
            "the filepath is './rnXXXXX.xml', where 'rnXXXXX' is derived from the specified map ID."
        )
    )
    args = parser.parse_args()
    return args

def get_map_id(args: argparse.Namespace) -> str:
    """
    Get the map ID from command line arguments.
    
    Parameters
    ==========
    args : argparse.Namespace
        Command line arguments.
    
    Returns
    =======
    str
        Map ID formatted 'mapXXXXX', where 'X' is a digit.
    """
    map_id = args.map_id
    assert re.match(map_id_pattern, map_id)
    return map_id

def get_map_filepath(args: argparse.Namespace) -> str:
    """
    Get the map download filepath from command line arguments.
    
    Parameters
    ==========
    args : argparse.Namespace
        Command line arguments.
        
    Returns
    =======
    str
        Download filepath for the 2x resolution map PNG file. Without the command line option
        specifying a filepath, the default is './mapXXXXX.png', where 'mapXXXXX' is the map ID.
    """
    map_filepath = args.map_filepath
    if map_filepath is None:
        map_filepath = os.path.join(os.getcwd(), f'{args.map_id}.png')
    else:
        assert os.access(os.path.dirname(map_filepath), os.W_OK)
    return map_filepath

def get_ko_kgml_filepath(args: argparse.Namespace) -> str:
    """
    Get the KO KGML filepath from command line arguments.
    
    Parameters
    ==========
    args : argparse.Namespace
        Command line arguments.
        
    Returns
    =======
    str
        Download filepath for the map KO KGML file. Without the command line option specifying a
        filepath, the default is './koXXXXX.xml', where 'koXXXXX' is derived from the map ID.
    """
    ko_kgml_filepath = args.ko_kgml_filepath
    if ko_kgml_filepath is None:
        ko_map_id = args.map_id.replace('map', 'ko')
        ko_kgml_filepath = os.path.join(os.getcwd(), f'{ko_map_id}.xml')
    else:
        assert os.access(os.path.dirname(ko_kgml_filepath), os.W_OK)
    return ko_kgml_filepath

def get_rn_kgml_filepath(args: argparse.Namespace) -> str:
    """
    Get the RN KGML filepath from command line arguments if setting up overview map files.
    
    Parameters
    ==========
    args : argparse.Namespace
        Command line arguments.
    
    Returns
    =======
    str
        Download filepath for the map RN KGML file. Without the command line option specifying a
        filepath, the default is './rnXXXXX.xml', where 'rnXXXXX' is derived from the map ID.
    """
    map_id = args.map_id
    is_overview_map = True if re.match(overview_map_id_pattern, map_id) else False
    if is_overview_map:
        rn_kgml_filepath = args.rn_kgml_filepath
        if rn_kgml_filepath is None:
            rn_map_id = args.map_id.replace('map', 'rn')
            rn_kgml_filepath = os.path.join(os.getcwd(), f'{rn_map_id}.xml')
        else:
            assert os.access(os.path.dirname(rn_kgml_filepath), os.W_OK)
    else:
        rn_kgml_filepath = None
    return rn_kgml_filepath

def download_map_image(map_id: str, map_filepath: str) -> None:
    """
    Download a KEGG pathway map as a 2x resolution PNG file.

    Parameters
    ==========
    map_id : str
        ID of map to download.
    
    map_filepath : str
        Download filepath for image.
    """
    print(f"Downloading {map_id} image to '{map_filepath}'")
    url = f'{kegg_rest_api_get}/{map_id}/image2x'
    urllib.request.urlretrieve(url, map_filepath)
    
def download_kgml(map_id: str, kgml_filepath: str, type: str) -> None:
    """
    Download a KGML file with elements of type 'ko' or 'rn' for the given pathway. This file is needed
    for map customization.
    
    Parameters
    ==========
    map_id : str
        Map ID.
    
    kgml_filepath : str
        Download KGML filepath.
    
    type : str
        Type of annotation in the KGML file, either 'ko' or 'rn'.
    """
    assert type in ('ko', 'rn')
    print(f"Downloading {map_id} KGML {type.upper()} file to '{kgml_filepath}'")
    kgml_id = map_id.replace('map', type)
    url = f'{kegg_rest_api_get}/{kgml_id}/kgml'
    urllib.request.urlretrieve(url, kgml_filepath)

def rescale_kgml(map_id: str, kgml_filepath: str, type: str, factor: float) -> None:
    """
    Rescale the sizes and positions of map objects (boxes, circles, lines) represented in a KGML
    file, overwriting the existing file.
    
    Parameters
    ==========
    map_id : str
        Map ID.
    
    kgml_filepath : str
        KGML filepath.
    
    type : str
        Type of annotation in the KGML file, either 'ko' or 'rn'.

    factor : float
        Rescaling factor.
    
    Returns
    =======
    None
    """
    assert type in ('ko', 'rn')
    assert isinstance(factor, float)
    assert factor > 0

    print(
        f"Rescaling {map_id} KGML {type.upper()} file by "
        f"{int(factor) if factor.is_integer() else factor}x"
    )
    
    with open(kgml_filepath) as f:
        kgml_text = f.read()

    params = ('x', 'y', 'width', 'height', 'coords')
    is_overview_map = True if re.match(overview_map_id_pattern, map_id) else False

    for param in params:
        if param == 'coords' and not is_overview_map:
            # Reaction lines are only found in overview maps and are parameterized by coordinates.
            continue
        
        # Go through the KGML text, looking for the parameter. Build a new KGML file in chunks
        # between each occurrence of the parameter.
        if param == 'coords':
            pattern = re.compile(r'coords="(.*)"')
            def get_new_string(groups: Tuple[str]) -> str:
                new_string = ''
                for value in groups[0].split(','):
                    new_value = float(value) * factor
                    new_value = int(new_value) if new_value.is_integer() else new_value
                    new_string += str(new_value) + ','
                new_string = new_string[: -1]
                return new_string
        else:
            pattern = re.compile(param + r'="(\d+)"')
            def get_new_string(groups : Tuple[str]) -> str:
                new_value = float(groups[0]) * factor
                new_value = int(new_value) if new_value.is_integer() else new_value
                new_string = str(new_value)
                return new_string
        
        new_kgml_text: str = None
        prev_chunk_stop = None
        for match in re.finditer(pattern, kgml_text):
            groups = match.groups()
            new_string = get_new_string(groups)
            if prev_chunk_stop:
                # Add the text from the previous pattern match through the rescaled match.
                new_kgml_text += kgml_text[prev_chunk_stop: match.start(1)] + new_string
            else:
                # This is the first pattern match.
                new_kgml_text = kgml_text[: match.start(1)] + new_string
            prev_chunk_stop = match.end(1)

        if prev_chunk_stop:
            # Add the text after the final pattern match.
            new_kgml_text += kgml_text[prev_chunk_stop: ]
        else:
            # There were no pattern matches.
            new_kgml_text = kgml_text

        kgml_text = new_kgml_text
    
    with open(kgml_filepath, 'w') as f:
        f.write(kgml_text)

def add_line_width(map_id: str, kgml_filepath: str, type: str, width: float = 12.0) -> None:
    """
    Add a width attribute to overview map line objects represented in a KGML file, overwriting the
    existing file.
    
    Parameters
    ==========
    map_id : str
        Overview map ID.
    
    kgml_filepath : str
        KGML filepath.
    
    type : str
        Type of annotation in the KGML file, either 'ko' or 'rn'.

    width : float, 12.0
        Line width. The default covers up lines in the underlying base map when rendering the KGML.
    
    Returns
    =======
    None
    """
    is_overview_map = True if re.match(overview_map_id_pattern, map_id) else False
    if not is_overview_map:
        return
    
    assert type in ('ko', 'rn')
    assert isinstance(width, float)
    assert width > 0

    print(f"Adding line width parameter to {map_id} KGML {type.upper()} file")
    
    with open(kgml_filepath) as f:
        kgml_text = f.read()
    
    width = int(width) if width.is_integer() else width
    kgml_text = kgml_text.replace('type="line" coords="', f'type="line" width="{width}" coords="')
    
    with open(kgml_filepath, 'w') as f:
        f.write(kgml_text)

if __name__ == '__main__':
    args = get_args()
    map_id = get_map_id(args)
    map_filepath = get_map_filepath(args)
    ko_kgml_filepath = get_ko_kgml_filepath(args)
    rn_kgml_filepath = get_rn_kgml_filepath(args)

    download_map_image(map_id, map_filepath)
    download_kgml(map_id, ko_kgml_filepath, 'ko')
    if rn_kgml_filepath is not None:
        download_kgml(map_id, rn_kgml_filepath, 'rn')

    rescale_kgml(map_id, ko_kgml_filepath, 'ko', 2.0)
    add_line_width(map_id, ko_kgml_filepath, 'ko')
    if rn_kgml_filepath is not None:
        rescale_kgml(map_id, rn_kgml_filepath, 'rn', 2.0)
        add_line_width(map_id, rn_kgml_filepath, 'rn')

    print(f"{map_id} files are set up")
