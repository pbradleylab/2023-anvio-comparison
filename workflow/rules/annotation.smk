include: "sample_selection.smk"

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

rule microbeannotator:
    input:
        profile=rules.untar_kofams_profile.output,
        ko_list=rules.untar_kofams_list.output,
        db=rules.microbeannotator_db_builder.output,
        ref=rules.transeq_cds.output,
    output: "results/annotation/microbeAnnotator/{genome}/kofam_results/{genome}.faa.kofam.filt"
    params:
        dir=directory("results/annotation/microbeAnnotator/{genome}/"),
        method="blast"
    conda: "../envs/microbeannotator.yml"
    log: "logs/annotation/microbeAnnotator/{genome}.log"
    shell:
        """
        mkdir -p {params.dir}
        microbeannotator -i {input.ref} -d {input.db} -o {params.dir} -m {params.method} 2> {log}
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

rule anvio_estimate_metabolism:
    input:
        done=rules.anvio_gen_contigs_db.output.done,
        kegg=rules.anvio_run_kegg_kofams.output,
        kofam=rules.anvio_setup_kegg_kofams.output
    output:"results/annotation/anvio/anvio_estimate_metabolism/{genome}/anvio_estimate_metabolism.0"
    params:
        db=rules.anvio_gen_contigs_db.output.db
    conda:"../envs/anvio.yml"
    log:"logs/annotation/anvio_estimate_metabolism/{genome}.log"
    shell:
        """
        anvi-estimate-metabolism -c {params.db} --kegg-data-dir {input.kofam} -o {output} 2> {log} 
        """

rule anvio_display_metabolism:
    input:
       done=rules.anvio_gen_contigs_db.output.done,
       metabol=rules.anvio_estimate_metabolism.output
    output:"/tmp/{genome}/anvio-display-metabolism.0"
    params:
       db=rules.anvio_gen_contigs_db.output.db
    conda:"../envs/anvio.yml"
    log:"logs/annotation/anvio-display-metabolism/{genome}.log"
    shell:
        """
        anvi-display-metabolism -c {params.db} 2> {log}
        """

rule create_enzyme_file_kofam:
    input:rules.kofamscan.output
    output:"results/annotation/create_enzyme_file/{genome}.kofam"
    shell:
        """
        echo "enzyme\tsource\torthology" > {output}
        cut -f1 {input} | sed 's/$/\tkofam/g' >> {output}
        """

rule create_enzyme_file_microbeannotator:
    input:rules.microbeannotator.output
    output:"results/annotation/create_enzyme_file/{genome}.microbeannotator"
    shell:
        """
        echo "enzyme\tsource\torthology" > {output}
        cut -f1 {input} | sed 's/$/\tKOfam/g' >> {output}
        """

# Useless for me
rule microbeannotator_script_gen_user_module_file:
    input:rules.create_enzyme_file_microbeannotator.output
    output:"results/annotation/gen_user_module_file/{genome}.microbeannotator"
    params:
       name="{genome}",
       categorization="User modules; kofam set; {genome} metabolism",
    conda:"../envs/anvio.yml"
    log:"logs/annotation/microbeannotator_script_gen_user_module_file/{genome}.log"
    shell:
        """
        list=$(cut -f1 {input} | paste -s -d+ -)
        anvi-script-gen-user-module-file -I {params.name} \
                  -n “Anvio comparison pathway analysis” \
                  -c {params.categorization} \
                  -e {input} \
                  -d $list” \
                  -o {output} 2> {log}
        """

rule kofam_script_gen_user_module_file:
    input:rules.create_enzyme_file_kofam.output
    output:"results/annotation/gen_user_module_file/{genome}.kofam"
    params:
       name="{genome}",
       categorization="User modules; kofam set; {genome} metabolism",
    conda:"../envs/anvio.yml"
    log:"logs/annotation/kofam_script_gen_user_module_file/{genome}.log"
    shell:
        """
        list=$(cut -f1 {input} | paste -s -d+ -)
        anvi-script-gen-user-module-file -I {params.name} \
                  -n “Anvio comparison pathway analysis” \
                  -c {params.categorization} \
                  -e {input} \
                  -d $list” \
                  -o {output} 2> {log}
        """
# Useless
rule microbeannotator_setup_user_modules:
    input:rules.microbeannotator_script_gen_user_module_file.output
    output:"results/annotation/microbeannotator_setup_user_modules/{genome}/"
    conda:"../envs/anvio.yml"
    log:"logs/annotation/microbeannotator_setup_user_modules/{genome}.log"
    shell:
        """
        anvi-setup-user-modules --user-modules {input} -o {output} 2> {log}
        """

rule kofam_setup_user_modules:
    input:rules.kofam_script_gen_user_module_file.output
    output:"results/annotation/kofam_setup_user_modules/{genome}/"
    conda:"../envs/anvio.yml"
    log:"logs/annotation/kofam_setup_user_modules/{genome}.log"
    shell:
        """
        anvi-setup-user-modules --user-modules {input} -o {output} 2> {log}
        """

rule microbeannotator_estimate_metabolism:
    input:rules.create_enzyme_file_microbeannotator.output
    output:"results/annotation/microbeannotator_setup_user_modules/{genome}/"
    params:
       db=rules.anvio_gen_contigs_db.output.db
    conda:"../envs/anvio.yml"
    log:"logs/annotation/microbeannotator_estimate_metabolism/{genome}.log"
    shell:
        """
        anvi-estimate-metabolism -c {params.db} --enzymes-txt {input} 2> {log} 
        """

rule kofam_estimate_metabolism:
    input:rules.kofam_setup_user_modules.output
    output:"results/annotation/kofam_setup_user_modules/{genome}/"
    params:
       db=rules.anvio_gen_contigs_db.output.db
    conda:"../envs/anvio.yml"
    log:"logs/annotation/kofam_estimate_metabolism/{genome}.log"
    shell:
        """
        anvi-estimate-metabolism -c {params.db} --user-modules {input} 2> {log}
        """

rule microbeannotator_display_metabolism:
    input:
    output:"results/annotation/microbeannotator_display_metabolism/{genome}_met.png"
    params:
       db=rules.anvio_gen_contigs_db.output.db
    conda:"../envs/anvio.yml"
    log:"logs/annotation/microbeannotator_display_metabolism/{genome}.log"
    shell:
        """
        anvi-display-metabolism -c {params.db} --file {output} 2> {log}
        """

rule kofam_display_metabolism:
    input:
    output:"results/annotation/kofam_display_metabolism/{genome}_met.png"
    params:
       db=rules.anvio_gen_contigs_db.output.db
    conda:"../envs/anvio.yml"
    log:"logs/annotation/kofam_display_metabolism/{genome}.log"
    shell:
        """
        anvi-display-metabolism -c {params.db} --file {output} 2> {log}
        """
