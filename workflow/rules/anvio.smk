include: "sample_selection.smk"
include: "setup.smk"


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
        db="results/annotation/anvio/contigs_db/anvio/default/{genome}/output.db",
        done="results/temporary/anvio/{genome}.anvio_gen_contigs_db"
    conda:"../envs/anvio.yml"
    log: "logs/annotation/contigs_db/anvio/default/{genome}.log"
    params:
        bacteria="{genome}"
    shell:
        """
        anvi-gen-contigs-database -f {input} -o {output.db} -n {params.bacteria} 2> {log}
        touch {output.done}
        """

rule anvio_make_gff:
    input:
        db=rules.anvio_gen_contigs_db.output.db
    output:"results/annotation/anvio/anvio_make_gff/{genome}.gff",
    conda:"../envs/anvio.yml"
    log: "logs/annotation/anvio_make_gff/{genome}.log"
    params:
        bacteria="{genome}"
    shell:
        """
        anvi-get-sequences-for-gene-calls -c {input.db} -o {output} --export-gff3
        """

rule anvio_gen_contigs_no_heuristic_db:
    input:rules.anvio_script_reformat.output
    output:
        db="results/annotation/anvio/contigs_db/anvio/no_hueristic/{genome}/output.db",
        done="results/temporary/anvio/{genome}.anvio_gen_contigs_no_heuristic_db"
    conda:"../envs/anvio.yml"
    log: "logs/annotation/contigs_db/anvio/default/{genome}.log"
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
    output:"results/temporary/anvio_run_kegg_kofams/{genome}/anvio_run_kegg_kofams.0"
    params:
       db=rules.anvio_gen_contigs_db.output.db
    conda:"../envs/anvio.yml"
    log: "logs/annotation/anvio_run_kegg_kofams/{genome}.log"
    threads: 40
    shell:
        """
        anvi-run-kegg-kofams -c {params.db} --kegg-data-dir /home/kananen.13/workflows/2023-anvio-comparison/resources/feb_anvio/no_stray -T {threads} --just-do-it 2> {log}
        touch {output}
        """

rule anvio_run_kegg_kofams_no_heuristic:
    input:
        done=rules.anvio_gen_contigs_no_heuristic_db.output.done,
        kofam=rules.anvio_setup_kegg_kofams.output
    output:"results/temporary/kegg/{genome}/anvio_run_kegg_kofams_no_heuristic/anvio_run_kegg_kofams.0"
    params:
       db=rules.anvio_gen_contigs_no_heuristic_db.output.db
    conda:"../envs/anvio.yml"
    log: "logs/annotation/anvio_run_kegg_kofams/no_heuristic/{genome}.log"
    threads: 40
    shell:
        """
        anvi-run-kegg-kofams --skip-bitscore-heuristic -c {params.db} --kegg-data-dir {input.kofam} -T {threads} --just-do-it 2> {log}
        touch {output}
        """

rule anvio_export_functions:
     input:
        done=rules.anvio_gen_contigs_db.output.done,
        kegg=rules.anvio_run_kegg_kofams.output
     output:"results/annotation/anvio/functions/default/{genome}.tsv"
     params:
        db=rules.anvio_gen_contigs_db.output.db
     conda:"../envs/anvio.yml"
     log: "logs/annotation/anvio/functions/default/{genome}.log"
     shell:
        """
        anvi-export-functions -c {params.db} -o {output} 2> {log}
        """

rule anvio_export_functions_no_huerestic:
     input:
        done=rules.anvio_gen_contigs_no_heuristic_db.output.done,
        kegg=rules.anvio_run_kegg_kofams_no_heuristic.output
     output:"results/annotation/anvio/functions/no_hueristic/{genome}.tsv"
     params:
        db=rules.anvio_gen_contigs_db.output.db
     conda:"../envs/anvio.yml"
     log: "logs/annotation/anvio/functions/no_hueristic/{genome}.log"
     shell:
        """
        anvi-export-functions -c {params.db} -o {output} 2> {log}
        """

rule anvio_gen_contigs_stray_db:
    input:rules.anvio_script_reformat.output
    output:
        db="results/annotation/anvio/contigs_db/anvio/stray/{genome}/output.db",
        done="results/temporary/anvio/{genome}.anvio_gen_contigs_stray_db"
    conda:"../envs/anvio_stray.yml"
    log: "logs/annotation/anvio/contigs_db/anvio/stray/{genome}.log"
    params:
        bacteria="{genome}"
    shell:
        """
        anvi-gen-contigs-database -f {input} -o {output.db} -n {params.bacteria} 2> {log}
        touch {output.done}
        """

rule anvio_run_kegg_kofams_stray:
    input:
        done=rules.anvio_gen_contigs_stray_db.output.done,
        kofam=rules.anvio_setup_kegg_stray_kofams.output
    output:"results/temporary/{genome}/anvio_run_kegg_kofams_stray/anvio_run_kegg_kofams.0"
    params:
       db=rules.anvio_gen_contigs_stray_db.output.db
    conda:"../envs/anvio_stray.yml"
    log: "logs/annotation/anvio_run_kegg_kofams/stray/{genome}.log"
    threads: 40
    shell:
        """
        anvi-run-kegg-kofams --include-stray-KOs -c {params.db} --kegg-data-dir {input.kofam} -T {threads} --just-do-it 2> {log}
        touch {output}
        """

rule anvio_export_functions_stray:
     input:
        done=rules.anvio_gen_contigs_stray_db.output.done,
        kegg=rules.anvio_run_kegg_kofams_stray.output
     output:"results/annotation/anvio/functions/stray/{genome}.tsv"
     params:
        db=rules.anvio_gen_contigs_stray_db.output.db
     conda:"../envs/anvio_stray.yml"
     log: "logs/annotation/anvio/functions/stray/{genome}.log"
     shell:
        """
        anvi-export-functions -c {params.db} -o {output} 2> {log}
        """

rule anvio_estimate_metabolism:
    input:
        done=rules.anvio_gen_contigs_db.output.done,
        kegg=rules.anvio_run_kegg_kofams.output,
        kofam=rules.anvio_setup_kegg_kofams.output
    output:"results/annotation/anvio/metabolism/default/{genome}/anvio_estimate_metabolism_modules.txt"
    params:
        db=rules.anvio_gen_contigs_db.output.db,
        prefix="results/annotation/anvio/metabolism/default/{genome}/anvio_estimate_metabolism"
    conda:"../envs/anvio.yml"
    log:"logs/annotation/anvio_estimate/metabolism/default/{genome}.log"
    shell:
        """
        anvi-estimate-metabolism -c {params.db} --kegg-data {input.kofam} -O {params.prefix} --just-do-it 2> {log} 
        """

rule anvio_estimate_metabolism_stray:
    input:
        done=rules.anvio_gen_contigs_stray_db.output.done,
        kegg=rules.anvio_run_kegg_kofams_stray.output,
        kofam=rules.anvio_setup_kegg_stray_kofams.output
    output:"results/annotation/anvio/metabolism/stray/{genome}/anvio_estimate_metabolism_modules.txt"
    params:
        db=rules.anvio_gen_contigs_stray_db.output.db,
        prefix="results/annotation/anvio/metabolism/stray/{genome}/anvio_estimate_metabolism"
    conda:"../envs/anvio_stray.yml"
    log:"logs/annotation/anvio/metabolism/stray/{genome}.log"
    shell:
        """
        anvi-estimate-metabolism -c {params.db} --kegg-data-dir {input.kofam} -O {params.prefix} --include-stray-KOs --just-do-it 2> {log}
        """

rule anvio_get_sequences_for_gene_calls:
    input:
        done=rules.anvio_gen_contigs_stray_db.output.done,
    output:"results/annotation/anvio/get_sequences_for_gene_calls/{genome}.faa"
    params:
        db=rules.anvio_gen_contigs_stray_db.output.db,
    conda:"../envs/anvio_stray.yml"
    log:"logs/annotation/anvio/get_sequences_for_gene_calls/{genome}.log"
    shell:
        """
        anvi-get-sequences-for-gene-calls -c {params.db} --get-aa-sequences -o {output} 2> {log}
        """
