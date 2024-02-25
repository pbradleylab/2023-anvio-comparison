#!env bash

pip install hmmer
pip install attrs

git clone git@github.com:cruizperez/MicrobeAnnotator.git
cd MicrobeAnnotator && git checkout 9bbc5f6
export PATH=$PATH:$PWD/MicrobeAnnotator/bin/

