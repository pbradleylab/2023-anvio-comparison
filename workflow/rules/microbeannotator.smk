include: "sample_selection.smk"
include: "setup.smk"

rule microbeannotator:
    input:
        profile=rules.untar_kofams_profile.output,
        ko_list=rules.untar_kofams_list.output,
        db=rules.microbeannotator_db_builder.output,
        ref=rules.anvi_get_sequences_for_gene_calls.output,
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

rule microbeannotator_refined:
    input:
        profile=rules.untar_kofams_profile.output,
        ko_list=rules.untar_kofams_list.output,
        db="/home/kananen.13/workflows/2023-anvio-comparison/resources/microbeannotator_full/",
        #db=rules.microbeannotator_db_builder.output,
        ref=rules.anvi_get_sequences_for_gene_calls.output,
    output: "results/annotation/microbeAnnotator_refined/{genome}/kofam_results/{genome}.faa.kofam.filt"
    params:
        dir=directory("results/annotation/microbeAnnotator_refined/{genome}/"),
        method="diamond"
    conda: "../envs/microbeannotator.yml"
    threads: 40
    log: "logs/annotation/microbeAnnotator/refined/{genome}.log"
    shell:
        """
        mkdir -p {params.dir}
        microbeannotator --threads {threads} --refine -i {input.ref} -d {input.db} -o {params.dir} -m {params.method} 2> {log}
        """

rule create_enzyme_file_microbeannotator:
    input:rules.microbeannotator.output
    output:"results/annotation/create_enzyme_file/{genome}.microbeannotator"
    shell:
        """
        echo "gene_id\tenzyme_accession\tsource" > {output}
        while read line
        do
            echo $line | cut -f1,3 -d' ' | sed 's/$/\tKOfam/g' | sed 's/ /\t/' >> /tmp/{wildcards.genome}.create_enzyme_file_microbeannotator
        done < {input}
        grep -v "#" /tmp/{wildcards.genome}.create_enzyme_file_microbeannotator >> {output}
        """

rule anvio_gen_contigs_no_annotations_refined_microbeannotator:
    input:rules.anvio_script_reformat.output
    output:
        db="results/annotation/anvio/anvio_gen_contigs_no_annotations_refined_microbeannotator/{genome}/output.db",
        done="/tmp/anvio/{genome}.anvio_gen_contigs_no_annotation_refined_microbeannotator"
    conda:"../envs/anvio.yml"
    log: "logs/annotation/anvio_gen_contigs_no_annotations_refined_microbeannotator/{genome}.log"
    params:
        bacteria="{genome}"
    shell:
        """
        anvi-gen-contigs-database -f {input} -o {output.db} -n {params.bacteria} 2> {log}
        touch {output.done}
        """

rule anvio_gen_contigs_no_annotations_microbeannotator:
    input:rules.anvio_script_reformat.output
    output:
        db="results/annotation/anvio/anvio_gen_contigs_no_annotations_microbeannotator/{genome}/output.db",
        done="/tmp/anvio/{genome}.anvio_gen_contigs_no_annotation_microbeannotator"
    conda:"../envs/anvio.yml"
    log: "logs/annotation/anvio_gen_contigs_no_annotations_microbeannotator/{genome}.log"
    params:
        bacteria="{genome}"
    shell:
        """
        anvi-gen-contigs-database -f {input} -o {output.db} -n {params.bacteria} 2> {log}
        touch {output.done}
        """

rule microbeannotator_estimate_metabolism:
    input:
       enzymes=rules.create_enzyme_file_microbeannotator.output,
       kofam=rules.anvio_setup_kegg_kofams.output,
       db=rules.anvio_gen_contigs_no_annotations_microbeannotator.output.db
    output:"results/annotation/microbeannotator_estimate_metabolism/{genome}/anvio_estimate_metabolism_modules.txt"
    params:
       prefix="results/annotation/microbeannotator_estimate_metabolism/{genome}/anvio_estimate_metabolism"
    conda:"../envs/anvio.yml"
    log:"logs/annotation/microbeannotator_estimate_metabolism/{genome}.log"
    shell:
        """
        anvi-estimate-metabolism --include-kos-not-in-kofam -c {input.db} --enzymes-txt {input.enzymes} --kegg-data-dir {input.kofam} -O {params.prefix} 2> {log}
        """

rule create_enzyme_file_microbeannotator_refined:
    input:rules.microbeannotator_refined.output
    output:"results/annotation/create_enzyme_file/{genome}.microbeannotator_refined"
    shell:
        """
        echo "gene_id\tenzyme_accession\tsource" > {output}
        while read line
        do
            echo $line | cut -f1,3 -d' ' | sed 's/$/\tKOfam/g' | sed 's/ /\t/' >> /tmp/{wildcards.genome}.create_enzyme_file_microbeannotator
        done < {input}
        grep -v "#" /tmp/{wildcards.genome}.create_enzyme_file_microbeannotator >> {output}
        """

rule microbeannotator_estimate_metabolism_refined:
    input:
       enzymes=rules.create_enzyme_file_microbeannotator_refined.output,
       kofam=rules.anvio_setup_kegg_kofams.output,
       db=rules.anvio_gen_contigs_no_annotations_refined_microbeannotator.output.db
    output:"results/annotation/microbeannotator_refined_estimate_metabolism/{genome}/anvio_estimate_metabolism_modules.txt"
    params:
       prefix="results/annotation/microbeannotator_refined_estimate_metabolism/{genome}/anvio_estimate_metabolism"
    conda:"../envs/anvio.yml"
    log:"logs/annotation/microbeannotator_refined_estimate_metabolism/{genome}.log"
    shell:
        """
        anvi-estimate-metabolism --include-kos-not-in-kofam -c {input.db} --enzymes-txt {input.enzymes} --kegg-data-dir {input.kofam} -O {params.prefix} 2> {log}
        """
