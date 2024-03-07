#!env bash

pip install -U git+https://github.com/cruizperez/MicrobeAnnotator@master
pip install hmmer
pip install attrs
pip install dataclasses
wget -c https://github.com/cruizperez/MicrobeAnnotator/blob/master/microbeannotator/data/01.KEGG_Regular_Module_Information.pickle
mv 01.KEGG_Regular_Module_Information.pickle $(echo $CONDA_PREFIX)/lib/python3.7/site-packages/microbeannotator/data/
