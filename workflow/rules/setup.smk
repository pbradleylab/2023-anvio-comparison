"""
Fetches the correct versions of the databases needed for running the analysis.

Author: Kathryn Kananen
"""

# Metadata file used for family classification
rule download_gtdb_meta:
   output: "resources/gtdb/bac120_taxonomy_r214.tsv"
   shell:
       """
       wget -c https://data.gtdb.ecogenomic.org/releases/release214/214.0/bac120_taxonomy_r214.tsv -O {output}
       """

# Databases needed for DRAM.
# Note dram wasn't utilized for this analysis 
rule download_dram_kegg_list:
    output: "resources/DRAM_data/kofam_ko_list.tsv"
    conda: "../envs/dram.yml"
    threads: 40
    shell:
        """
        DRAM-setup.py prepare_databases --threads 40 --select_db 'kofam_ko_list' --output resources/DRAM_data
        """
rule download_dram_hmm:
    output:"resources/DRAM_data/kofam_profiles.hmm"
    conda: "../envs/dram.yml"
    threads: 40
    shell: "../envs/dram.yml"
        """
        DRAM-setup.py prepare_databases --threads {threads} --select_db 'kofam_hmm' --output {output}
        """

# KOfam database and ko_list
rule download_kofams_profile:
    output: "resources/kofams/profile.tar.gz"
    shell:
        """
        wget -c --no-http-keep-alive ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz -O {output}
        """
rule untar_kofams_profile:
    input: rules.download_kofams_profile.output
    output: directory("resources/kofams/profiles")
    shell:
        """
        tar -xvzf {input}
        mv ./profiles $(dirname {output})
        """
rule download_kofams_list:
    output: "resources/kofams/ko_list.gz"
    shell:
        """
        wget -c --no-http-keep-alive ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz -O {output}
        """
rule untar_kofams_list:
    input: rules.download_kofams_list.output
    output: "resources/kofams/list/ko_list"
    shell:
        """
        gunzip -c {input} > {output}
        """

# Download the eggnog databases used for checking strays
rule eggnog_db:
    output:
        outdir=directory("resources/eggnog/db/"),
        db="resources/eggnog/db/eggnog.db"
    conda: "../envs/comparison.yml"
    shell:
        """
        download_eggnog_data.py --data_dir {output.outdir} -y
        """
rule eggnog_bacteria_db:
    output:
        dmnd="resources/eggnog/bacteria/db/bacteria.dmnd",
        outdir=directory("resources/eggnog/bacteria/db/")
    conda: "../envs/comparison.yml"
    shell:
        """
        mkdir -p {output.outdir}
        create_dbs.py -y -m diamond --dbname bacteria --taxa Bacteria --data_dir {output.outdir}
        """

# Create two separate databases for each version of anvio, one without and one with strays.
rule anvio_setup_kegg_kofams:
    output:directory("resources/anvio/no_stray/")
    conda:"../envs/anvio_stray.yml"
    threads: 40
    shell:
        """
        anvi-setup-kegg-data -T {threads} --download-from-kegg --kegg-data-dir {output}
        """
rule anvio_setup_kegg_stray_kofams:
    output:directory("resources/anvio/stray/")
    conda:"../envs/anvio_stray.yml"
    threads: 40
    shell:
        """
        anvi-setup-kegg-data -T {threads} --download-from-kegg --include-stray-KOs --kegg-data-dir {output}
        """

# Build the microbeannotator database.
rule microbeannotator_db_builder:
    output: directory("resources/microbeannotator/")
    conda: "../envs/microbeannotator.yml"
    params:
        method="diamond"
    shell:
        """
        microbeannotator_db_builder --database {output} -m {params.method} --no_aspera
        """
