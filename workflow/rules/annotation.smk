include: "sample_selection.smk"


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
    output:directory("resources/anvio/kofams/")
    conda:"../envs/anvio.yml"
    shell:
        """
        anvi-setup-kegg-data --kegg-data-dir {output}
        touch {output}
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


rule kofamscan:
    input:
        ko_list=rules.untar_kofams_list.output,
        profiles=rules.untar_kofams_profile.output,
        faa=rules.transeq.output
    output:"results/annotation/kofamscan/{genome}.tsv"
    conda:"../envs/kofamscan.yml"
    shell:
        """
        exec_annotation -f detail-tsv -p {input.profiles} -o {output} {input.faa} -k {input.ko_list} --tmp-dir /tmp/{wildcards.genome}
        """

rule kofamscan_params:
    input:
        ko_list=rules.untar_kofams_list.output,
        profiles=rules.untar_kofams_profile.output,
        faa=rules.transeq.output
    output:"results/annotation/kofamscan_params_run/{genome}.tsv"
    conda:"../envs/kofamscan.yml"
    shell:
        """
        exec_annotation -f detail-tsv -p {input.profiles} -o {output} {input.faa} -k {input.ko_list} --tmp-dir /tmp/params_run/{wildcards.genome}
        """

rule microbeannotator:
    input:
        profile=rules.untar_kofams_profile.output,
        ko_list=rules.untar_kofams_list.output,
        db=rules.microbeannotator_db_builder.output,
        ref=rules.transeq_cds.output,
    output: "results/annotation/microbeAnnotator/{genome}/kofam_results/{genome}_cds.faa.kofam.filt"
    params:
        dir=directory("results/annotation/microbeAnnotator/{genome}/"),
        method="blast"
    conda: "../envs/microbeannotator.yml"
    log: "logs/annotation/microbeAnnotator/{genome}.log"
    shell:
        """
        mkdir -p {params.dir} >> microbeannotator_comms.sh
        ls {output} "||" /home/kananen.13/workflows/tools/MicrobeAnnotator/bin/microbeannotator -i {input.ref} -d {input.db} -o {params.dir} -m {params.method} >> microbeannotator_comms.sh 2> {log}
        """

rule microbeannotator_refined:
    input:
        profile=rules.untar_kofams_profile.output,
        ko_list=rules.untar_kofams_list.output,
        db=rules.microbeannotator_db_builder.output,
        ref=rules.transeq_cds.output,
    output: "results/annotation/microbeAnnotator/refined/{genome}/kofam_results/{genome}.faa.kofam.filt"
    params:
        dir=directory("results/annotation/microbeAnnotator/refined/{genome}/"),
        method="blast"
    conda: "../envs/microbeannotator.yml"
    log: "logs/annotation/microbeAnnotator/refined/{genome}.log"
    shell:
        """
        mkdir -p {params.dir} >> make_dir.sh
        #echo microbeannotator --refine -i {input.ref} -d {input.db} -o {params.dir} -m {params.method} >> microbeannotator_comms.sh 2> {log}
        """

rule anvio_script_reformat:
    input: rules.download_genomes.output
    output:"results/annotation/anvio/reformat/{genome}.fasta"
    conda:"../envs/anvio.yml"
    params:
    log: "logs/annotation/anvio_script_reformat_fasta/{genome}.log"
    shell:
        """
        anvi-script-reformat-fasta {input} -o {output} --simplify-name --seq-type NT 2> {log}
        """

rule anvio_gen_contigs_db:
    input:rules.anvio_script_reformat.output
    output:
        db="results/annotation/anvio/anvio_gen_contigs_db/{genome}/output.db",
        done="/tmp/anvio/{genome}.anvio_gen_contigs_db"
    conda:"../envs/anvio.yml"
    log: "logs/annotation/anvio_gen_contigs_db/{genome}.log"
    params:
        bacteria="{genome}"
    shell:
        """
        anvi-gen-contigs-database -f {input} -o {output.db} -n {params.bacteria} 2> {log}
        touch {output.done}
        """

rule anvio_gen_contigs_no_heuristic_db:
    input:rules.anvio_script_reformat.output
    output:
        db="results/annotation/anvio/anvio_gen_contigs_no_heuristic_db/{genome}/output.db",
        done="/tmp/anvio/{genome}.anvio_gen_contigs_no_heuristic_db"
    conda:"../envs/anvio.yml"
    log: "logs/annotation/anvio_gen_contigs_db/{genome}.log"
    params:
        bacteria="{genome}"
    shell:
        """
        anvi-gen-contigs-database -f {input} -o {output.db} -n {params.bacteria} 2> {log}
        touch {output.done}
        """

rule anvio_run_kegg_kofams:
    input:
        done=rules.anvio_gen_contigs_db.output.done,
        kofam=rules.anvio_setup_kegg_kofams.output
    output:"/tmp/{genome}/anvio_run_kegg_kofams.0"
    params:
       db=rules.anvio_gen_contigs_db.output.db
    conda:"../envs/anvio.yml"
    log: "logs/annotation/anvio_run_kegg_kofams/{genome}.log"
    threads: 1
    shell:
        """
        anvi-run-kegg-kofams -c {params.db} --kegg-data-dir {input.kofam} -T {threads} --just-do-it 2> {log}
        touch {output}
        """

rule anvio_run_kegg_kofams_no_heuristic:
    input:
        done=rules.anvio_gen_contigs_no_heuristic_db.output.done,
        kofam=rules.anvio_setup_kegg_kofams.output
    output:"/tmp/{genome}/anvio_run_kegg_kofams_no_heuristic/anvio_run_kegg_kofams.0"
    params:
       db=rules.anvio_gen_contigs_db.output.db
    conda:"../envs/anvio.yml"
    log: "logs/annotation/anvio_run_kegg_kofams/no_heuristic/{genome}.log"
    threads: 1
    shell:
        """
        anvi-run-kegg-kofams --skip-bitscore-heuristic -c {params.db} --kegg-data-dir {input.kofam} -T {threads} --just-do-it 2> {log}
        touch {output}
        """

rule anvio_export_functions:
     input:
        done=rules.anvio_gen_contigs_db.output.done,
        kegg=rules.anvio_run_kegg_kofams.output
     output:"results/annotation/anvio/anvio_functions/{genome}.tsv"
     params:
        db=rules.anvio_gen_contigs_db.output.db
     conda:"../envs/anvio.yml"
     log: "logs/annotation/anvio_export_functions/{genome}.log"
     shell:
        """
        anvi-export-functions -c {params.db} -o {output} 2> {log}
        """

rule anvio_export_functions_no_huerestic:
     input:
        done=rules.anvio_gen_contigs_no_heuristic_db.output.done,
        kegg=rules.anvio_run_kegg_kofams_no_heuristic.output
     output:"results/annotation/anvio/anvio_functions_no_huerestic/{genome}.tsv"
     params:
        db=rules.anvio_gen_contigs_db.output.db
     conda:"../envs/anvio.yml"
     log: "logs/annotation/anvio_export_functions_no_huerestic/{genome}.log"
     shell:
        """
        anvi-export-functions -c {params.db} -o {output} 2> {log}
        """

rule dram:
    input:
        hmm=rules.download_dram_hmm.output,
        ko=rules.download_dram_kegg_list.output,
        ref=rules.anvio_script_reformat.output
    output:
        gff="results/annotation/dram/{genome}/genes.gff",
        annotations="results/annotation/dram/{genome}/annotations.tsv",
        rrna="results/annotation/dram/{genome}/rrnas.tsv",
        trna="results/annotation/dram/{genome}/trnas.tsv"
    threads: 40
    params:
        outdir=directory("results/annotation/dram/{genome}")
    log: "logs/annotation/dram/{genome}.log"
    conda: "../envs/dram.yml"
    shell:
        """
        rm -r /tmp/dram/ || echo "making dram folder /tmp/dram/"
        DRAM-setup.py import_config --config_loc config/dram.og.json
        DRAM.py annotate -i {input.ref} -o /tmp/dram/ --threads {threads} 2> {log}
        mkdir -p {params.outdir} && mv /tmp/dram/* {params.outdir}
        """

rule anvio_estimate_metabolism:
    input:
        done=rules.anvio_gen_contigs_db.output.done,
        kegg=rules.anvio_run_kegg_kofams.output,
        kofam=rules.anvio_setup_kegg_kofams.output
    output:"results/annotation/anvio/anvio_estimate_metabolism/{genome}/anvio_estimate_metabolism_modules.txt"
    params:
        db=rules.anvio_gen_contigs_db.output.db,
        prefix="results/annotation/anvio/anvio_estimate_metabolism/{genome}/anvio_estimate_metabolism"
    conda:"../envs/anvio.yml"
    log:"logs/annotation/anvio_estimate_metabolism/{genome}.log"
    shell:
        """
        anvi-estimate-metabolism -c {params.db} --kegg-data-dir {input.kofam} -O {params.prefix} 2> {log} 
        """

rule create_enzyme_file_kofam:
    input:rules.kofamscan.output
    output:"results/annotation/create_enzyme_file/{genome}.kofam"
    shell:
        """
        echo "gene_id\tenzyme_accession\tsource" > {output}
        while read line
        do
            echo $line | cut -f1,3 -d' ' | grep -v '#' | sed 's/$/\tKOfam/g' | sed 's/ /\t/' >> {output}
        done < {input}
        """

rule create_enzyme_file_microbeannotator:
    input:rules.microbeannotator.output
    output:"results/annotation/create_enzyme_file/{genome}.microbeannotator"
    shell:
        """
        echo "gene_id\tenzyme_accession\tsource > {output}
        while read line
        do
            echo $line | cut -f1,3 -d' ' | grep -v '#' | sed 's/$/\tKOfam/g' | sed 's/ /\t/' >> {output}
        done < {input}
        """

rule microbeannotator_estimate_metabolism:
    input:
       enzymes=rules.create_enzyme_file_microbeannotator.output,
       kofam=rules.anvio_setup_kegg_kofams.output
    output:"results/annotation/microbeannotator_estimate_metabolism/{genome}/anvio_estimate_metabolism_modules.txt"
    params:
       db=rules.anvio_gen_contigs_db.output.db,
       prefix="results/annotation/microbeannotator_estimate_metabolism/{genome}/anvio_estimate_metabolism"
    conda:"../envs/anvio.yml"
    log:"logs/annotation/microbeannotator_estimate_metabolism/{genome}.log"
    shell:
        """
        anvi-estimate-metabolism --include-kos-not-in-kofam -c {params.db} --enzymes-txt {input.enzymes} --kegg-data-dir {input.kofam} -O {params.prefix} 2> {log} 
        """

rule kofam_estimate_metabolism:
    input:
       enzymes=rules.create_enzyme_file_kofam.output,
       kofam=rules.anvio_setup_kegg_kofams.output
    output: "results/annotation/kofam_estimate_metabolism/{genome}/anvio_estimate_metabolism_modules.txt"
    params:
       db=rules.anvio_gen_contigs_db.output.db,
       prefix="results/annotation/microbeannotator_estimate_metabolism/{genome}/anvio_estimate_metabolism"
    conda:"../envs/anvio.yml"
    log:"logs/annotation/kofam_estimate_metabolism/{genome}.log"
    shell:
        """
        anvi-estimate-metabolism -c {params.db} --user-modules {input.enzymes} --kegg-data-dir {input.kofam} -O {params.prefix} 2> {log}
        """

