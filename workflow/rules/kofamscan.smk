include: "sample_selection.smk"
include: "setup.smk"
import glob


rule kofamscan:
    input:
        ko_list=rules.untar_kofams_list.output,
        profiles=rules.untar_kofams_profile.output,
        indir=rules.anvi_get_sequences_for_gene_calls.output
    output:"results/annotation/kofamscan/{genome}.tsv"
    conda:"../envs/kofamscan.yml"
    threads:40
    shell:
        """
        exec_annotation --cpu {threads} -f detail-tsv -p {input.profiles} -o {output} {input.indir} -k {input.ko_list} --tmp-dir /tmp/{wildcards.genome}/kofamscan/
        """

rule kofamscan_refined:
    input:
        ko_list=rules.untar_kofams_list.output,
        profiles=rules.untar_kofams_profile.output,
        faa=rules.anvi_get_sequences_for_gene_calls.output
    output:"results/annotation/kofamscan_params_run/{genome}.tsv"
    conda:"../envs/kofamscan.yml"
    shell:
        """
        exec_annotation -T 0.5 -E 0.00001 -f detail-tsv -p {input.profiles} -o {output} {input.faa} -k {input.ko_list} --tmp-dir /tmp/params_run/{wildcards.genome}
        """

rule create_enzyme_file_kofam:
    input:rules.kofamscan.output
    output:"results/annotation/create_enzyme_file/{genome}.kofam"
    shell:
        """
        echo "gene_id\tenzyme_accession\tsource" > {output}
        grep '*' {input} > /tmp/{wildcards.genome}.tsv
        cut -f2,3 /tmp/{wildcards.genome}.tsv > /tmp/{wildcards.genome}.1
        awk '{{print $0"\tKOfam"}}' /tmp/{wildcards.genome}.1 >> {output}
        """

rule anvio_gen_contigs_no_annotations_refined_kofam:
    input:rules.anvio_script_reformat.output
    output:
        db="results/annotation/anvio/anvio_gen_contigs_no_annotations_refined_kofam/{genome}/output.db",
        done="/tmp/anvio/{genome}.anvio_gen_contigs_no_annotation_refined"
    conda:"../envs/anvio.yml"
    log: "logs/annotation/anvio_gen_contigs_no_annotations_refined_kofam/{genome}.log"
    params:
        bacteria="{genome}"
    shell:
        """
        anvi-gen-contigs-database -f {input} -o {output.db} -n {params.bacteria} 2> {log}
        touch {output.done}
        """

rule anvio_gen_contigs_no_annotations_kofam:
    input:rules.anvio_script_reformat.output
    output:
        db="results/annotation/anvio/anvio_gen_contigs_no_annotations_kofam/{genome}/output.db",
        done="/tmp/anvio/{genome}.anvio_gen_contigs_no_annotation_kofam"
    conda:"../envs/anvio.yml"
    log: "logs/annotation/anvio_gen_contigs_no_annotations_kofam/{genome}.log"
    params:
        bacteria="{genome}"
    shell:
        """
        anvi-gen-contigs-database -f {input} -o {output.db} -n {params.bacteria} 2> {log}
        touch {output.done}
        """

rule kofam_estimate_metabolism:
    input:
       enzymes=rules.create_enzyme_file_kofam.output,
       kofam=rules.anvio_setup_kegg_kofams.output,
       db=rules.anvio_gen_contigs_no_annotations_kofam.output.db
    output: "results/annotation/kofam_estimate_metabolism/{genome}/anvio_estimate_metabolism_modules.txt"
    params:
       prefix="results/annotation/kofam_estimate_metabolism/{genome}/anvio_estimate_metabolism"
    conda:"../envs/anvio.yml"
    log:"logs/annotation/kofam_estimate_metabolism/{genome}.log"
    shell:
        """
        anvi-estimate-metabolism -c {input.db} --enzymes-txt {input.enzymes} --kegg-data-dir {input.kofam} -O {params.prefix} 2> {log}
        """
