sed '1,3d' microbeannotator/Bradyrhizobium_manausense_BR3351/kofam_results/Bradyrhizobium_manausense_BR3351.faa.kofam.filt | cut -f3 | sort -u > comparison/microbeannotator/Bradyrhizobium_manausense_BR3351.kofam
sed '1,3d' microbeannotator/Lactiplantibacillus_plantarum_SRCM100442/kofam_results/Lactiplantibacillus_plantarum_SRCM100442.faa.kofam.filt | cut -f3 | sort -u > comparison/microbeannotator/Lactiplantibacillus_plantarum_SRCM100442.kofam
sed '1,3d' microbeannotator/Robinsoniella_peoriensis_B23985/kofam_results/Robinsoniella_peoriensis_B23985.faa.kofam.filt | cut -f3 | sort -u > comparison/microbeannotator/Robinsoniella_peoriensis_B23985.kofam
sed '1,3d' microbeannotator/Vibrio_cholerae_RFB16/kofam_results/Vibrio_cholerae_RFB16.faa.kofam.filt | cut -f3 | sort -u > comparison/microbeannotator/Vibrio_cholerae_RFB16.kofam

sed '1d' anvio/functions/Bradyrhizobium_manausense_BR3351.txt | grep KOfam | cut -f3 | sort -u > comparison/anvio/Bradyrhizobium_manausense_BR3351.txt
sed '1d' anvio/functions/Lactiplantibacillus_plantarum_SRCM100442.txt | grep KOfam | cut -f3 | sort -u > comparison/anvio/Lactiplantibacillus_plantarum_SRCM100442.txt
sed '1d' anvio/functions/Robinsoniella_peoriensis_B23985.txt | grep KOfam | cut -f3 | sort -u > comparison/anvio/Robinsoniella_peoriensis_B23985.txt
sed '1d' anvio/functions/Vibrio_cholerae_RFB16.txt | grep KOfam | cut -f3 | sort -u > comparison/anvio/Vibrio_cholerae_RFB16.txt

grep -F -v -x -f comparison/microbeannotator/Bradyrhizobium_manausense_BR3351.kofam comparison/anvio/Bradyrhizobium_manausense_BR3351.txt > comparison/in_anvio_only/Bradyrhizobium_manausense_BR3351.txt
grep -F -v -x -f comparison/microbeannotator/Lactiplantibacillus_plantarum_SRCM100442.kofam comparison/anvio/Lactiplantibacillus_plantarum_SRCM100442.txt > comparison/in_anvio_only/Lactiplantibacillus_plantarum_SRCM100442.txt
grep -F -v -x -f comparison/microbeannotator/Robinsoniella_peoriensis_B23985.kofam comparison/anvio/Robinsoniella_peoriensis_B23985.txt > comparison/in_anvio_only/Robinsoniella_peoriensis_B23985.txt
grep -F -v -x -f comparison/microbeannotator/Vibrio_cholerae_RFB16.kofam comparison/anvio/Vibrio_cholerae_RFB16.txt > comparison/in_anvio_only/Vibrio_cholerae_RFB16.txt

grep -F -v -x -f comparison/anvio/Bradyrhizobium_manausense_BR3351.txt comparison/microbeannotator/Bradyrhizobium_manausense_BR3351.kofam > comparison/in_microbeannotator_only/Bradyrhizobium_manausense_BR3351.txt
grep -F -v -x -f comparison/anvio/Lactiplantibacillus_plantarum_SRCM100442.txt comparison/microbeannotator/Lactiplantibacillus_plantarum_SRCM100442.kofam > comparison/in_microbeannotator_only/Lactiplantibacillus_plantarum_SRCM100442.txt
grep -F -v -x -f comparison/anvio/Robinsoniella_peoriensis_B23985.txt comparison/microbeannotator/Robinsoniella_peoriensis_B23985.kofam > comparison/in_microbeannotator_only/Robinsoniella_peoriensis_B23985.txt
grep -F -v -x -f comparison/anvio/Vibrio_cholerae_RFB16.txt comparison/microbeannotator/Vibrio_cholerae_RFB16.kofam > comparison/in_microbeannotator_only/Vibrio_cholerae_RFB16.txt
