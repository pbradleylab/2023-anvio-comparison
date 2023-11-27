# Assumes that you are in an enviroment where anvio v8 is installed
anvi-setup-kegg-data

anvi-script-reformat-fasta --simplify-names --output ./anvio/reformat/Bradyrhizobium_manausense_BR3351.fna genomes/Bradyrhizobium_manausense_BR3351.fna
anvi-script-reformat-fasta --simplify-names --output ./anvio/reformat/Robinsoniella_peoriensis_B23985.fna genomes/R_peoriensis_B23985_assembly.fna
anvi-script-reformat-fasta --simplify-names --output ./anvio/reformat/Lactiplantibacillus_plantarum_SRCM100442.fna genomes/Lactiplantibacillus_plantarum_SRCM100442.fna
anvi-script-reformat-fasta --simplify-names --output ./anvio/reformat/Vibrio_cholerae_RFB16.fna genomes/Vibrio_cholerae_RFB16.fna

anvi-gen-contigs-database -f anvio/reformat/Bradyrhizobium_manausense_BR3351.fna -o anvio/dbs/Bradyrhizobium_manausense_BR3351.db
anvi-gen-contigs-database -f anvio/reformat/Robinsoniella_peoriensis_B23985.fna -o anvio/dbs/Robinsoniella_peoriensis_B23985.db
anvi-gen-contigs-database -f anvio/reformat/Lactiplantibacillus_plantarum_SRCM100442.fna -o anvio/dbs/Lactiplantibacillus_plantarum_SRCM100442.db
anvi-gen-contigs-database -f anvio/reformat/Vibrio_cholerae_RFB16.fna -o anvio/dbs/Vibrio_cholerae_RFB16.db

anvi-run-kegg-kofams -c anvio/dbs/Bradyrhizobium_manausense_BR3351.db
anvi-run-kegg-kofams -c anvio/dbs/Robinsoniella_peoriensis_B23985.db
anvi-run-kegg-kofams -c anvio/dbs/Lactiplantibacillus_plantarum_SRCM100442.db
anvi-run-kegg-kofams -c anvio/dbs/Vibrio_cholerae_RFB16.db

  
anvi-export-functions -c anvio/dbs/Bradyrhizobium_manausense_BR3351.db -o anvio/functions/Bradyrhizobium_manausense_BR3351.txt
anvi-export-functions -c anvio/dbs/Robinsoniella_peoriensis_B23985.db  -o anvio/functions/Robinsoniella_peoriensis_B23985.txt
anvi-export-functions -c anvio/dbs/Lactiplantibacillus_plantarum_SRCM100442.db  -o anvio/functions/Lactiplantibacillus_plantarum_SRCM100442.txt
anvi-export-functions -c anvio/dbs/Vibrio_cholerae_RFB16.db  -o anvio/functions/Vibrio_cholerae_RFB16.txt
