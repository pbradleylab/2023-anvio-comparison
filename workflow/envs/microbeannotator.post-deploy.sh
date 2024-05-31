#!env bash

pip install -U git+https://github.com/cruizperez/MicrobeAnnotator@master
pip install hmmer
pip install attrs
pip install dataclasses

wget -c https://github.com/cruizperez/MicrobeAnnotator/raw/master/microbeannotator/data/01.KEGG_Regular_Module_Information.pickle --no-check-certificate
wget -c https://github.com/cruizperez/MicrobeAnnotator/raw/master/microbeannotator/data/03.KEGG_Structural_Module_Information.pickle --no-check-certificate
wget -c https://github.com/cruizperez/MicrobeAnnotator/raw/master/microbeannotator/data/02.KEGG_Bifurcating_Module_Information.pickle --no-check-certificate

wget -c https://raw.githubusercontent.com/cruizperez/MicrobeAnnotator/master/microbeannotator/data/01.KEGG_DB/00.KEGG_Data_Scrapper.py --no-check-certificate
wget -c https://raw.githubusercontent.com/cruizperez/MicrobeAnnotator/master/microbeannotator/data/01.KEGG_DB/00.Module_Names.txt --no-check-certificate
wget -c https://raw.githubusercontent.com/cruizperez/MicrobeAnnotator/master/microbeannotator/data/01.KEGG_DB/01.Bifurcating_List.txt --no-check-certificate
wget -c https://raw.githubusercontent.com/cruizperez/MicrobeAnnotator/master/microbeannotator/data/01.KEGG_DB/02.Structural_List.txt --no-check-certificate
wget -c https://raw.githubusercontent.com/cruizperez/MicrobeAnnotator/master/microbeannotator/data/01.KEGG_DB/03.Bifurcating_Modules.dict --no-check-certificate
wget -c https://raw.githubusercontent.com/cruizperez/MicrobeAnnotator/master/microbeannotator/data/01.KEGG_DB/04.Structural_Modules.dict --no-check-certificate
wget -c https://raw.githubusercontent.com/cruizperez/MicrobeAnnotator/master/microbeannotator/data/01.KEGG_DB/05.Modules_Parsed.txt --no-check-certificate
wget -c https://raw.githubusercontent.com/cruizperez/MicrobeAnnotator/master/microbeannotator/data/01.KEGG_DB/06.Module_Groups.txt --no-check-certificate

mv 01.KEGG_Regular_Module_Information.pickle $(echo $CONDA_PREFIX)/lib/python3.7/site-packages/microbeannotator/data/
mv 02.KEGG_Bifurcating_Module_Information.pickle $(echo $CONDA_PREFIX)/lib/python3.7/site-packages/microbeannotator/data/
mv 03.KEGG_Structural_Module_Information.pickle $(echo $CONDA_PREFIX)/lib/python3.7/site-packages/microbeannotator/data/

mkdir $(echo $CONDA_PREFIX)/lib/python3.7/site-packages/microbeannotator/data/01.KEGG_DB/
mv 00.KEGG_Data_Scrapper.py $(echo $CONDA_PREFIX)/lib/python3.7/site-packages/microbeannotator/data/01.KEGG_DB/
mv 00.Module_Names.txt $(echo $CONDA_PREFIX)/lib/python3.7/site-packages/microbeannotator/data/01.KEGG_DB/
mv 01.Bifurcating_List.txt $(echo $CONDA_PREFIX)/lib/python3.7/site-packages/microbeannotator/data/01.KEGG_DB/
mv 02.Structural_List.txt $(echo $CONDA_PREFIX)/lib/python3.7/site-packages/microbeannotator/data/01.KEGG_DB/
mv 03.Bifurcating_Modules.dict $(echo $CONDA_PREFIX)/lib/python3.7/site-packages/microbeannotator/data/01.KEGG_DB/
mv 04.Structural_Modules.dict $(echo $CONDA_PREFIX)/lib/python3.7/site-packages/microbeannotator/data/01.KEGG_DB/
mv 05.Modules_Parsed.txt $(echo $CONDA_PREFIX)/lib/python3.7/site-packages/microbeannotator/data/01.KEGG_DB/
mv 06.Module_Groups.txt $(echo $CONDA_PREFIX)/lib/python3.7/site-packages/microbeannotator/data/01.KEGG_DB/
