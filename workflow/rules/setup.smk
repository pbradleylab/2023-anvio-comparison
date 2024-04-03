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

rule download_kofams_profile:
    output: "resources/kofams/profile.tar.gz"
    shell:
        """
        wget -c --no-http-keep-alive ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz -O {output}
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

rule untar_kofams_profile:
    input: rules.download_kofams_profile.output
    output: directory("resources/kofams/profiles")
    shell:
        """
        tar -xvzf {input}
        mv ./profiles $(dirname {output})
        """

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

rule microbeannotator_db_builder:
    output: directory("resources/microbeannotator/")
    conda: "../envs/microbeannotator.yml"
    params:
        method="diamond"
    shell:
        """
        microbeannotator_db_builder --database {output} -m {params.method} --no_aspera
        """
