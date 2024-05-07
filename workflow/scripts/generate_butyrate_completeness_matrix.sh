#!/bin/bash
# a script to estimate pathwise completeness of a custom butyrate biosynthesis module in Lachnospiraceae genomes 
#   using input annotation sets from each annotation tool
# designed to work with the output of the Snakemake workflow
# creates a matrix of the results wherein each row describes a genome and each column describes an annotation tool


# create user-defined modules database with the butyrate biosynthesis module
# modules file must be stored at this relative path: `Butanoate_module/modules/BUTANOATE``
anvi-setup-user-modules -u Butanoate_module/

# convert functions output from each tool to enzymes-txt input files for anvi-estimate-metabolism
# relies on the data being stored at these relative paths: `Archive/{TOOL}/functions/default/{GENOME_ACCESSION}.tsv`
# also relies on the input file `Lachno_genomes.txt` containing the accessions of the Lachnospiraceae genomes
mkdir -p Archive/enzymes-txt-files

## first we process the anvi'o files
mkdir -p Archive/enzymes-txt-files/anvio
mkdir -p Archive/enzymes-txt-files/anvio/default
while read acc; do \
  python functions_output_to_enzymes_txt.py Archive/anvio/functions/default/${acc}.tsv \
    Archive/enzymes-txt-files/anvio/default/${acc}_enzymes.txt \
    anvio; \
done < <(tail -n+2 Lachno_genomes.txt | cut -f 11 )

## then kofamscan files
mkdir -p Archive/enzymes-txt-files/kofamscan
mkdir -p Archive/enzymes-txt-files/kofamscan/default
while read acc; do \
  python functions_output_to_enzymes_txt.py Archive/kofamscan/functions/default/${acc}.tsv \
    Archive/enzymes-txt-files/kofamscan/default/${acc}_enzymes.txt \
    kofamscan; \
done < <(tail -n+2 Lachno_genomes.txt | cut -f 11 )

## then MicrobeAnnotator files
mkdir -p Archive/enzymes-txt-files/MicrobeAnnotator
mkdir -p Archive/enzymes-txt-files/MicrobeAnnotator/default
while read acc; do \
  python functions_output_to_enzymes_txt.py Archive/MicrobeAnnotator/functions/default/${acc}/annotation_results/${acc}.faa.annot \
    Archive/enzymes-txt-files/MicrobeAnnotator/default/${acc}_enzymes.txt \
    MicrobeAnnotator; \
done < <(tail -n+2 Lachno_genomes.txt | cut -f 11 )

# estimate completeness of the pathway
## first for anvi'o
mkdir -p Archive/BUTANOATE_MODULE_OUTPUT
mkdir -p Archive/BUTANOATE_MODULE_OUTPUT/anvio
mkdir -p Archive/BUTANOATE_MODULE_OUTPUT/anvio/default
while read acc; do \
  anvi-estimate-metabolism --enzymes-txt Archive/enzymes-txt-files/anvio/default/${acc}_enzymes.txt \
    -u Butanoate_module/ \
   --only-user-modules \
   -O Archive/BUTANOATE_MODULE_OUTPUT/anvio/default/${acc} \
   --include-zeros; \
done < <(tail -n+2 Lachno_genomes.txt | cut -f 11 )

## then for kofamscan
mkdir -p Archive/BUTANOATE_MODULE_OUTPUT/kofamscan
mkdir -p Archive/BUTANOATE_MODULE_OUTPUT/kofamscan/default
while read acc; do \
  anvi-estimate-metabolism --enzymes-txt Archive/enzymes-txt-files/kofamscan/default/${acc}_enzymes.txt \
    -u Butanoate_module/ \
   --only-user-modules \
   -O Archive/BUTANOATE_MODULE_OUTPUT/kofamscan/default/${acc} \
   --include-zeros; \
done < <(tail -n+2 Lachno_genomes.txt | cut -f 11 )

## then for microbeannotator
mkdir -p Archive/BUTANOATE_MODULE_OUTPUT/MicrobeAnnotator
mkdir -p Archive/BUTANOATE_MODULE_OUTPUT/MicrobeAnnotator/default
while read acc; do \
  anvi-estimate-metabolism --enzymes-txt Archive/enzymes-txt-files/MicrobeAnnotator/default/${acc}_enzymes.txt \
    -u Butanoate_module/ \
   --only-user-modules \
   -O Archive/BUTANOATE_MODULE_OUTPUT/MicrobeAnnotator/default/${acc} \
   --include-zeros; \
done < <(tail -n+2 Lachno_genomes.txt | cut -f 11 )

# finally, combine the results into one matrix
python gen_butyrate_matrix.py
