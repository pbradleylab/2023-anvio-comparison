include: "sample_selection.smk"
include: "setup.smk"

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
        mkdir -p {params.dir}
        microbeannotator -i {input.ref} -d {input.db} -o {params.dir} -m {params.method} 2> {log}
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
        mkdir -p {params.dir}
        microbeannotator --refine -i {input.ref} -d {input.db} -o {params.dir} -m {params.method} 2> {log}
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
            grep -v "#" >> {output}
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
        anvi-estimate-metabolism --include-kos-not-in-kofam -c {params.db} --enzymes-txt {input.enzymes} --kegg-data-dir /home/kananen.13/workflows/2023-anvio-comparison/resources/feb_anvio/stray -O {params.prefix} 2> {log}
        """
