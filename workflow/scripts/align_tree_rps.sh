#!/bin/bash

# first argument is directory containing fasta files, optional second argument is Cdd database

if [ $2 -eq 0 ]; then
  CDD_DIR=resources/cdd/Cdd
else 
  CDD_DIR=$2
fi

if [ -d $1 ]; then
  cd $1
else
  echo >&2 "Error: $1 does not exist in the current path or is not a directory"

cd $1

# align
mkdir -p aln
cd aln
for x in ../*faa; do echo $x; clustalo -i $x -o `basename $x .faa`.aln; done
cd ..

# build trees
mkdir -p tree
cd tree
for x in ../*faa; do echo $x; fasttree -out `basename $x .faa`.tree $x; done
cd ..

# run rpsblast
mkdir -p cdd
cd cdd
for i in `ls ../*.faa`; do echo $i
  rpsblast -db $CDD_DIR -query $i -evalue 1e-4 -out $(basename $i .faa).txt -outfmt 6
done
cd ..
