#!/usr/bin/python3
"""
Draw a customized KEGG map, selecting KOs to color by category in pathway maps, and reactions and
compounds associated with KOs to color in overview maps.

Author: Samuel Miller
Contact: samuelmiller10@gmail.com
"""

import os
import re
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

from io import StringIO
from itertools import combinations
from Bio.KEGG.KGML import KGML_parser
from typing import Dict, List, Set, Tuple
from Bio.Graphics.KGML_vis import KGMLCanvas

map_id_pattern = re.compile(r'map\d{5}')
# Overview maps have a certain range of IDs.
overview_map_id_pattern = re.compile(r'map\d{1}1[12]\d{2}')

ortholog_element_name_pattern = re.compile(r'ko:(K\d{5})')

# Uncategorized elements in overview maps are colored light gray.
overview_map_uncategorized_color = '#EDEDED'
# The diameter of compound circles is scaled up to cover up the underlying circles in the base map.
# The factor is based on the default KGML diameter in the KGML file of 28.
overview_map_circle_scaling_factor = 1.18

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
        help="ID of the map to color, formatted as 'mapXXXXX', where 'X' is a digit."
    )
    parser.add_argument(
        'ko_categories',
        help=(
            "Tab-delimited file of KO IDs in the first column headed 'KO', and an arbitrary "
            "number of additional columns for each category, each headed with a unique name that "
            "cannot contain commas. Values in the table are True or False, indicating whether the "
            "category contains the KO."
        )
    )
    parser.add_argument(
        '--output_filepath',
        default=None,
        help=(
            "Path of the colored map PDF file. Without this option and without the map name "
            "option, the default construction of the path is '<working_dir>/<map_id>.pdf'. With a "
            "map name, the default construction is '<working_dir>/<map_id>-<map_name>.pdf'."
        )
    )
    parser.add_argument(
        '--map_name',
        default=None,
        help="Name to give to the map, used in the default construction of the output filepath."
    )
    parser.add_argument(
        '--map_filepath',
        default=None,
        help=(
            "Filepath of the map PNG file. By default, without this option, the filepath is "
            "'./mapXXXXX.png', where 'mapXXXXX' is the specified map ID."
        )
    )
    parser.add_argument(
        '--ko_kgml_filepath',
        default=None,
        help=(
            "Filepath of the KGML KO file. By default, without this option, the filepath is "
            "'./koXXXXX.xml', where 'koXXXXX' is derived from the specified map ID."
        )
    )
    parser.add_argument(
        '--rn_kgml_filepath',
        default=None,
        help=(
            "Filepath of a KGML RN file, which is only needed for overview maps. By default, "
            "without this option, the filepath is './rnXXXXX.xml', where 'rnXXXXX' is derived from "
            "the specified map ID."
        )
    )
    parser.add_argument(
        '--without_base',
        action='store_true',
        default=False,
        help=(
            "Only draw elements defined by the KGML file and do not draw the base reference map "
            "image under these elements. This is especially useful for overview maps, as the "
            "reference map can include a small number of reactions and compounds not represented "
            "in KGML files and therefore not removable from the image, and these items and labels "
            "can have colors that clash with the new color scheme."
        )
    )
    parser.add_argument(
        '--category_colors',
        default=None,
        help=(
            "This option explicitly controls how categories and combinations of categories in the "
            "KO categories table are colored. Without this option, colors are automatically "
            "assigned to KOs using all categories and combinations of categories in the KO "
            "categories table and a default color palette. This option takes a tab-delimited file "
            "of two columns, the first column headed 'categories' and the second headed 'color'. "
            "Each value in the first column is a category with the same name as a column in the KO "
            "categories table, OR a combination of these categories, for KOs in multiple "
            "categories. Combinations of categories should be separated by a comma, e.g., "
            "categories A, B, and C, would be named 'A,B,C'. Each value in the second column is a "
            "color hex code. The order of rows of the table determines which color takes precedent "
            "on the map, with categories closer to the top of the table superceding categories "
            "closer to the bottom. Therefore, it makes sense for 'A,B,C' to be above 'A,B', 'A,C', "
            "and 'B,C', and for these to be above 'A', 'B', and 'C'. Lastly, note that "
            "combinations of categories must be explicit. If 'A,B' and 'A,C' are in the table but "
            "not 'B,C', then KOs in categories A plus B and A plus C are colored, but KOs in "
            "categories B plus C are not."
        )
    )
    parser.add_argument(
        '--only_color_shared',
        action='store_true',
        default=False,
        help=(
            "If True, and if colors are not explicitly given by the optional category colors "
            "table, then only color KOs that are shared by all categories in the KO categories "
            "table. If False, and if colors are not explicitly given by the optional category "
            "colors table, then KOs in any category or combination of categories are colored. "
            "Combinations of a greater number of categories are colored over fewer categories, "
            "e.g., if there are two categories, 'A' and 'B', then KOs in both categories will "
            "always receive the color for both rather than the color for either 'A' or 'B'."
        )
    )
    parser.add_argument(
        '--colormap',
        default='viridis',
        help=(
            "Name of a matplotlib colormap to use if colors are not explicitly given by the "
            "optional category colors table."
        )
    )
    parser.add_argument('--fontsize', type=int, default=18, help="Fontsize of map element labels.")
    
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

def load_ko_categories_table(args: argparse.Namespace) -> pd.DataFrame:
    """
    Load the KO categories table given command line arguments.
    
    Parameters
    ==========
    args : argparse.Namespace
        Command line arguments.
        
    Returns
    =======
    pd.DataFrame
        Table of categories of KO IDs.
    """
    ko_categories_filepath = args.ko_categories
    ko_categories_table = pd.read_csv(ko_categories_filepath, sep='\t', index_col=0)
    ko_categories_table.astype(bool)
    assert ko_categories_table.index.name == 'KO'
    assert len(ko_categories_table.columns)
    assert len(set(ko_categories_table.columns)) == len(ko_categories_table.columns)
    assert ko_categories_table.dtypes.eq(bool).all()
    for col in ko_categories_table.columns:
        assert ',' not in col
    return ko_categories_table

def get_map_filepath(args: argparse.Namespace) -> str:
    """
    Get the map filepath from command line arguments.
    
    Parameters
    ==========
    args : argparse.Namespace
        Command line arguments.
        
    Returns
    =======
    str
        Map filepath.
    """
    map_filepath = args.map_filepath
    if map_filepath is None:
        if args.without_base:
            return
        else:
            map_filepath = os.path.join(os.getcwd(), f'{args.map_id}.png')
    assert os.path.exists(map_filepath)
    return map_filepath

def get_ko_kgml_filepath(args: argparse.Namespace) -> str:
    """
    Get the path of a KO-annotated KGML file for the map from command line arguments.
    
    Parameters
    ==========
    args : argparse.Namespace
        Command line arguments.
        
    Returns
    =======
    str
        KO-annotated KGML filepath.
    """
    ko_kgml_filepath = args.ko_kgml_filepath
    if ko_kgml_filepath is None:
        ko_map_id = args.map_id.replace('map', 'ko')
        ko_kgml_filepath = os.path.join(os.getcwd(), f'{ko_map_id}.xml')
    assert os.path.exists(ko_kgml_filepath)
    return ko_kgml_filepath

def get_rn_kgml_filepath(args: argparse.Namespace) -> str:
    """
    Get the path of a RN-annotated KGML file for the overview map from command line arguments.
    
    Parameters
    ==========
    args : argparse.Namespace
        Command line arguments.

    Returns
    =======
    str
        RN-annotated KGML filepath.
    """
    map_id = args.map_id
    is_overview_map = True if re.match(overview_map_id_pattern, map_id) else False
    if not is_overview_map:
        return None
    rn_kgml_filepath = args.rn_kgml_filepath
    if rn_kgml_filepath is None:
        rn_map_id = map_id.replace('map', 'rn')
        rn_kgml_filepath = os.path.join(os.getcwd(), f'{rn_map_id}.xml')
    assert os.path.exists(rn_kgml_filepath)
    return rn_kgml_filepath

def get_output_filepath(args: argparse.Namespace, map_id: str, map_name: str = None) -> str:
    """
    Get the path of the colored map PDF output.
    
    Without an output path provided in the command line arguments, a default pathname is constructed
    as follows. Without a map name provided in the command line arguments, the construction is
    '<working_dir>/<map_id>.pdf'. With a map name provided in the command line arguments, the
    construction is '<working_dir>/<map_id>-<map_name>.pdf'.
    
    Parameters
    ==========
    args : argparse.Namespace
        Command line arguments.
    
    map_id : str
        The ID of the map of interest.
    
    map_name : str, None
        Map name to be incorporated into output filename.
    
    Returns
    =======
    str
        Output filepath.
    """
    output_filepath = args.output_filepath
    if output_filepath is None:
        if map_name is None:
            output_filepath = f"{map_id}.pdf"
        else:
            output_filepath = f"{map_id}-{map_name}.pdf"
    if os.path.isabs(output_filepath):
        assert os.access(os.path.dirname(output_filepath), os.W_OK)
    else:
        assert os.access(os.getcwd(), os.W_OK)
    return output_filepath

def get_category_colors_table(
    args: argparse.Namespace,
    ko_categories_table: pd.DataFrame
) -> pd.DataFrame:
    """
    Get a table of KO colors to be displayed given command line arguments.
    
    The table is either loaded from a file specified in the command line arguments or generated
    automatically.
    
    Parameters
    ==========
    args : argparse.Namespace
        Command line arguments.
    
    ko_categories_table : pd.DataFrame
        Table of categories of KO IDs.
        
    Returns
    =======
    pd.DataFrame
        Table of colors per category or combination of categories.
    """
    category_colors_filepath = args.category_colors
    colormap_name = args.colormap
    only_color_shared = args.only_color_shared

    if category_colors_filepath is not None:
        assert os.path.exists(category_colors_filepath)
        print(f"Loading color scheme from the category colors file, '{category_colors_filepath}'")
        category_colors_table = load_category_colors_table(
            category_colors_filepath, ko_categories_table
        )
        return category_colors_table
    
    # In the absence of a provided file, automatically generate a category colors table.
    categories = ko_categories_table.columns.to_list()
    if only_color_shared:
        category_combos: List[Tuple[str]] = list(combinations(categories, len(categories)))
    else:
        category_combos: List[Tuple[str]] = []
        # Larger combinations of categories are towards the top of the table and individual
        # categories are at the bottom to prioritize coloring larger combinations of categories.
        for combo_size in range(len(categories), 0, -1):
            category_combos += list(combinations(categories, combo_size))
    
    # Retrieve hex codes for each color drawn from the colormap.
    assert colormap_name in plt.colormaps()
    colormap = plt.cm.get_cmap(colormap_name, len(category_combos))
    hex_codes: List[str] = []
    for i in range(len(category_combos)):
        hex_codes.append(plt.colors.rgb2hex(colormap(i)))
    
    # The entries in the categories column are tuples of category names.
    pd.DataFrame.from_dict({'categories': category_combos, 'color': hex_codes})
    return category_colors_table
        
def load_category_colors_table(
    category_colors_filepath: str,
    ko_categories_table: pd.DataFrame
) -> pd.DataFrame:
    """
    Load a table of KO colors to be displayed.
    
    Parameters
    ==========
    category_colors_filepath : str
        Filepath to table of color assignments for categories and combinations of categories.
    
    ko_categories_table : pd.DataFrame
        Table of categories of KO IDs.
    
    Returns
    =======
    pd.DataFrame
        Table of colors per category or combination of categories.
    """
    category_colors_table = pd.read_csv(category_colors_filepath, sep='\t')
    assert category_colors_table.columns.to_list() == ['categories', 'color']
    assert len(set(category_colors_table['categories'])) == len(category_colors_table)
    
    possible_categories = ko_categories_table.columns.to_list()
    color_category_combos: List[Tuple[str]] = []
    for entry in category_colors_table['categories']:
        for category in entry.split(','):
            assert category in possible_categories
        color_category_combos.append(tuple(entry.split(',')))
    
    # The entries in the categories column are tuples of category names.
    category_colors_table = pd.DataFrame.from_dict(
        {'categories': color_category_combos, 'color': category_colors_table['color']}
    )
    return category_colors_table

def get_ko_colors(
    category_colors_table: pd.DataFrame,
    ko_categories_table: pd.DataFrame
) -> Dict[str, Tuple[str, int]]:
    """
    Map the ID of each KO to be colored to a hex code.
    
    The color is assigned by the highest priority category combination containing the KO.
    
    Parameters
    ==========
    category_colors_table : pd.DataFrame
        Table of colors per category (length 1 tuple) or combination of categories (longer tuple).
    
    ko_categories_table : pd.DataFrame
        Table of categories of KO IDs.

    Returns
    =======
    Dict[str, Tuple[str, int]]
        Keys are KO IDs and values are tuples of the color hex code and priority of the category
        combination represented by positive integers 1, 2, etc., with larger integers having higher
        priority.
    """
    ko_colors: Dict[str, Tuple[str, int]] = {}
    priority = len(category_colors_table)
    for category_combo, hex_code in zip(
        category_colors_table['categories'], category_colors_table['color']
    ):
        ko_category_combo_table = ko_categories_table[list(category_combo)]
        ko_ids: List[str] = ko_category_combo_table[
            ko_category_combo_table.all(axis=1)
        ].index.to_list()
        for ko_id in ko_ids:
            if ko_id in ko_colors:
                # The KO was already assigned a color in a higher priority category combination.
                continue
            ko_colors[ko_id] = (hex_code, priority)
        priority -= 1
    return ko_colors

def draw_map(
    map_id: str,
    ko_kgml_filepath: str,
    ko_colors: Dict[str, Tuple[str, int]],
    output_filepath: str,
    rn_kgml_filepath: str,
    map_filepath: str = None,
    map_name: str = None,
    without_base: bool = False,
    fontsize: int = 18
) -> None:
    """
    Draw a map with KOs colored by category.
    
    Parameters
    ==========
    map_id : str
        The ID of the map of interest.
    
    ko_kgml_filepath : str
        Path to a KO-annotated KGML file for the map.
    
    ko_colors : Dict[str, Tuple[str, int]]
        Map the ID of each KO to be colored to a hex code and priority of the category combination.
    
    output_filepath : str
        Path of the colored map PDF file.
    
    rn_kgml_filepath : str
        Path to a RN-annotated KGML file for an overview map.
    
    map_filepath : str, None
        Map filepath. This can only be None if the base reference map is not being drawn.
    
    map_name : str, None
        Map name to be incorporated into output filename.

    without_base : bool, False
        Draw the output image without the underlying reference map if True, with the map if False.
    
    fontsize : int, 18
        Fontsize of element labels from KGML files.
    
    Returns
    =======
    None
    """
    print(
        f"Drawing {map_id}{': ' + map_name if map_name else ''}"
        f"{' without base map' if without_base else ''}"
    )

    is_overview_map = True if re.match(overview_map_id_pattern, map_id) else False
    if is_overview_map:
        assert rn_kgml_filepath is not None
        draw_overview_map(
            ko_kgml_filepath,
            rn_kgml_filepath,
            ko_colors,
            output_filepath,
            map_filepath=map_filepath,
            without_base=without_base,
            fontsize=fontsize
        )
        return
    
    tree = ET.ElementTree(file=ko_kgml_filepath)
    root = tree.getroot()

    for element in root:
        if element.tag != 'entry':
            continue
        if element.attrib['type'] != 'ortholog':
            continue

        graphics_elements = [
            subelement for subelement in element if subelement.tag == 'graphics'
        ]
        assert len(graphics_elements) > 0
        
        element_name = element.attrib['name']
        matches = re.finditer(ortholog_element_name_pattern, element_name)

        prev_priority = 0
        priority = None
        for match in matches:
            ortholog_ko_id = match.group(1)
            try:
                hex_code, priority = ko_colors[ortholog_ko_id]
            except KeyError:
                continue
            
            if priority < prev_priority:
                continue
            prev_priority = priority

            for graphics_element in graphics_elements:
                graphics_element.attrib['bgcolor'] = hex_code

        if priority is None:
            for graphics_element in graphics_elements:
                # Uncolored KOs have a white background.
                graphics_element.attrib['bgcolor'] = '#FFFFFF'
    
    pathway = KGML_parser.read(
        StringIO(ET.tostring(root, encoding='unicode', method='xml'))
    )
    pathway.image = map_filepath
    
    canvas = KGMLCanvas(
        pathway,
        import_imagemap=not without_base,
        label_compounds=False,
        fontsize=fontsize
    )
    canvas.draw(output_filepath)

def draw_overview_map(
    ko_kgml_filepath: str,
    rn_kgml_filepath: str,
    ko_colors: Dict[str, Tuple[str, int]],
    output_filepath: str,
    map_filepath: str = None,
    without_base: bool = False,
    fontsize: int = 18
) -> None:
    """
    Draw an overview map with KOs colored by category.
    
    Overview maps have a certain range of KEGG IDs.

    Parameters
    ==========
    
    ko_kgml_filepath : str
        Path to a KO-annotated KGML file for the map.
        
    rn_kgml_filepath : str
        Path to a RN-annotated KGML file for the map.
    
    ko_colors : Dict[str, Tuple[str, int]]
        Map the ID of each KO to be colored to a hex code and priority of the category combination.
    
    output_filepath : str
        Path of the colored map PDF file.

    map_filepath : str, None
        Map filepath. This can only be None if the base reference map is not being drawn.

    without_base : bool, False
        Draw the output image without the underlying reference map if True, with the map if False.

    fontsize : int, 18
        Fontsize of element labels from KGML files.

    Returns
    =======
    None
    """
    ko_tree = ET.ElementTree(file=ko_kgml_filepath)
    ko_root = ko_tree.getroot()
    
    # Parse each ortholog entry (line), which corresponds to a reaction element. In order to
    # subsequently set the colors of compound circles associated with reactions, record the color
    # hex code and priority of the category combination containing the KO.
    element_order_dict: Dict[int, List[Tuple[str, ET.Element]]] = {}
    ortholog_entry_ids: List[str] = []
    reaction_colors: Dict[str, Tuple[str, int]] = {}
    for element in ko_root:
        if element.tag != 'entry':
            continue
        if element.attrib['type'] != 'ortholog':
            continue
        
        entry_id = element.attrib['id']
        ortholog_entry_ids.append(entry_id)
        
        graphics_elements = [
            subelement for subelement in element if subelement.tag == 'graphics'
        ]
        assert len(graphics_elements) > 0
        
        element_name = element.attrib['name']
        matches = re.finditer(ortholog_element_name_pattern, element_name)

        prev_priority = 0
        priority = None
        for match in matches:
            ortholog_ko_id = match.group(1)
            try:
                hex_code, priority = ko_colors[ortholog_ko_id]
            except KeyError:
                continue
            
            if priority < prev_priority:
                continue
            prev_priority = priority

            for graphics_element in graphics_elements:
                graphics_element.attrib['fgcolor'] = hex_code
        
        if priority is None:
            for graphics_element in graphics_elements:
                graphics_element.attrib['fgcolor'] = overview_map_uncategorized_color
            try:
                element_order_dict[-1].append((entry_id, element))
            except KeyError:
                element_order_dict[-1] = [(entry_id, element)]
            continue
        
        # Record the ortholog information for the corresponding reaction element.
        assert entry_id not in reaction_colors
        reaction_colors[entry_id] = (hex_code, priority)
        try:
            element_order_dict[priority].append((entry_id, element))
        except KeyError:
            element_order_dict[priority] = [(entry_id, element)]
    
    # Parse each reaction element, which references at least one substrate and at least one product.
    # For each compound entry (circle), record the color hex code and priority of the highest
    # priority category combination involving it.
    compound_entry_colors: Dict[str, Tuple[str, int]] = {}
    for element in ko_root:
        if element.tag != 'reaction':
            continue
        
        element_id = element.attrib['id']
        try:
            hex_code, priority = reaction_colors[element_id]
        except KeyError:
            continue
        
        substrate_elements = [
            subelement for subelement in element if subelement.tag == 'substrate'
        ]
        assert len(substrate_elements) > 0
        product_elements = [
            subelement for subelement in element if subelement.tag == 'product'
        ]
        assert len(product_elements) > 0

        for subelement in substrate_elements + product_elements:
            subelement_id = subelement.attrib['id']
            try:
                prev_value = compound_entry_colors[subelement_id]
            except KeyError:
                compound_entry_colors[subelement_id] = (hex_code, priority)
                continue
            
            prev_priority = prev_value[1]
            if priority > prev_priority:
                compound_entry_colors[subelement_id] = (hex_code, priority)
    
    # Parse each compound entry (circle). If the entry participates in colored KO reactions, set the
    # circle color to that of the highest priority category combination.
    ko_compound_entry_ids: List[str] = []
    for element in ko_root:
        if element.tag != 'entry':
            continue
        if element.attrib['type'] != 'compound':
            continue
        
        graphics_elements = [
            subelement for subelement in element if subelement.tag == 'graphics'
        ]
        assert len(graphics_elements) > 0
        
        element_id = element.attrib['id']
        ko_compound_entry_ids.append(element_id)
        try:
            hex_code, priority = compound_entry_colors[element_id]
            try:
                element_order_dict[priority].append((entry_id, element))
            except KeyError:
                element_order_dict[priority] = [(entry_id, element)]
        except KeyError:
            hex_code = overview_map_uncategorized_color
            try:
                element_order_dict[-1].append((entry_id, element))
            except KeyError:
                element_order_dict[-1] = [(entry_id, element)]
        
        for graphics_element in graphics_elements:
            graphics_element.attrib['fgcolor'] = hex_code
            graphics_element.attrib['bgcolor'] = hex_code
            graphics_element.attrib['width'] = str(
                float(graphics_element.attrib['width']) * overview_map_circle_scaling_factor
            )
            graphics_element.attrib['height'] = str(
                float(graphics_element.attrib['height']) * overview_map_circle_scaling_factor
            )
        
    ko_compound_entry_ids: Set[str] = set(ko_compound_entry_ids)

    # Fill in uncategorized reactions and compounds from the base map that are not associated with a
    # KO, and are therefore not represented in the KGML KO file, but are associated with a KEGG
    # REACTION, and are therefore represented in the KGML RN file. Each overview map has a KGML RN
    # file, unlike many maps for individual pathways.
    rn_tree = ET.ElementTree(file=rn_kgml_filepath)
    rn_root = rn_tree.getroot()

    reaction_entry_ids: List[str] = []
    for element in rn_root:
        if element.tag != 'entry':
            continue
        if element.attrib['type'] != 'reaction':
            continue

        entry_id = element.attrib['id']
        if entry_id in ortholog_entry_ids:
            continue

        graphics_elements = [
            subelement for subelement in element if subelement.tag == 'graphics'
        ]
        assert len(graphics_elements) > 0
        
        for graphics_element in graphics_elements:
            graphics_element.attrib['fgcolor'] = overview_map_uncategorized_color
        
        try:
            element_order_dict[-1].append((entry_id, element))
        except KeyError:
            element_order_dict[-1] = [(entry_id, element)]
        
        # Add the reaction element from the KGML RN XML to the tree made from the KGML KO XML.
        ko_root.append(element)
        reaction_entry_ids.append(entry_id)
    
    rn_compound_entry_ids: List[str] = []
    for element in rn_root:
        if element.tag != 'reaction':
            continue
        
        entry_id = element.attrib['id']
        if entry_id not in reaction_entry_ids:
            continue
        
        # The drawer requires that both the reaction entry element (already added) and this reaction
        # definition element are in the KGML.
        ko_root.append(element)
        
        substrate_elements = [
            subelement for subelement in element if subelement.tag == 'substrate'
        ]
        assert len(substrate_elements) > 0
        product_elements = [
            subelement for subelement in element if subelement.tag == 'product'
        ]
        assert len(product_elements) > 0
        
        for subelement in substrate_elements + product_elements:
            subelement_id = subelement.attrib['id']
            rn_compound_entry_ids.append(subelement_id)
        
        element_order_dict[-1].append((entry_id, element))
    rn_compound_entry_ids = set(rn_compound_entry_ids).difference(ko_compound_entry_ids)
    
    for element in rn_root:
        if element.tag != 'entry':
            continue
        if element.attrib['type'] != 'compound':
            continue
        
        entry_id = element.attrib['id']
        if entry_id not in rn_compound_entry_ids:
            continue
        
        graphics_elements = [
            subelement for subelement in element if subelement.tag == 'graphics'
        ]
        assert len(graphics_elements) > 0
        
        for graphics_element in graphics_elements:
            graphics_element.attrib['fgcolor'] = overview_map_uncategorized_color
            graphics_element.attrib['bgcolor'] = overview_map_uncategorized_color
            graphics_element.attrib['width'] = str(
                float(graphics_element.attrib['width']) * overview_map_circle_scaling_factor
            )
            graphics_element.attrib['height'] = str(
                float(graphics_element.attrib['height']) * overview_map_circle_scaling_factor
            )
        
        try:
            element_order_dict[-1].append((entry_id, element))
        except KeyError:
            element_order_dict[-1] = [(entry_id, element)]
        
        ko_root.append(element)

    new_root = ET.Element(ko_root.tag)
    for priority in sorted(element_order_dict):
        values = element_order_dict[priority]
        values = sorted(values, key=lambda value: value[0])
        for value in values:
            element = value[1]
            new_root.append(element)

    pathway = KGML_parser.read(
        StringIO(ET.tostring(new_root, encoding='unicode', method='xml'))
    )
    pathway.image = map_filepath
    
    canvas = KGMLCanvas(
        pathway,
        import_imagemap=not without_base,
        label_orthologs=False,
        label_reaction_entries=False,
        label_compounds=False,
        fontsize=fontsize
    )
    canvas.non_reactant_transparency = 1
    canvas.draw(output_filepath)

if __name__ == '__main__':
    args = get_args()
    map_id = get_map_id(args)
    ko_categories_table = load_ko_categories_table(args)
    map_name = args.map_name
    without_base = args.without_base
    map_filepath = get_map_filepath(args)
    ko_kgml_filepath = get_ko_kgml_filepath(args)
    rn_kgml_filepath = get_rn_kgml_filepath(args)
    output_filepath = get_output_filepath(args, map_id, map_name=map_name)
    category_colors_table = get_category_colors_table(args, ko_categories_table)
    fontsize = args.fontsize
    
    ko_colors = get_ko_colors(category_colors_table, ko_categories_table)
    category_combos: List[Tuple[str]] = category_colors_table['categories'].to_list()
    draw_map(
        map_id,
        ko_kgml_filepath,
        ko_colors,
        output_filepath,
        map_filepath=map_filepath,
        rn_kgml_filepath=rn_kgml_filepath,
        map_name=map_name,
        without_base=without_base,
        fontsize=fontsize
    )
