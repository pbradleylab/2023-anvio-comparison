#include:"sample_selection.smk"
from os import walk


def get_protein_fasta(wildcards):
    out = []
    dir = rules.symlink_genomes.output
    genomes = next(walk(dir), (None, None, []))[2]
    proteins = [x for x in genomes if x.endswith("protein.faa.gz")]
    for protein in proteins:
        out.append(rules.anvio_export_functions.output[1].format(sample=protein))
    return(out)

def get_protein_fasta(wildcards):
    out = []
    dir = rules.symlink_genomes.output
    genomes = next(walk(dir), (None, None, []))[2]
    nucleotides = [x for x in genomes if x.endswith("genomic.fna.gz")]
    for nucleotide in nucleotides:
        out.append(rules.anvio_export_functions.output[1].format(sample=protien))
    return(out)


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
        gzip -dv {input} > {output}
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
        microbeannotator_db_builder --light --database {output} -m {params.method} --no_aspera
        """

rule kofamscan:
    input:
        ko_list=rules.download_kofams_list.output,
        profiles=rules.untar_kofams_profile.output,
        faa=rules.symlink_genomes.output
    output:"results/annotation/kofamscan/{genome}.tsv"
    shell:
        """
        ./exec_annotation -f mapper -p {input.faa} -o {output} {input.faa} -k {input.ko_list}
        """




rule microbeannotator:
    input:
        profile=rules.untar_kofams_profile.output,
        ko_list=rules.untar_kofams_list.output,
        db=rules.microbeannotator_db_builder.output,
        ref=get_protein_fasta
    output: "results/{project}/annotation/microbeAnnotator/{subsample}/kofam_results/{subsample}.faa.kofam.filt"
    params:
        dir=directory("results/{project}/annotation/microbeAnnotator/{subsample}/"),
        method="blast"
    conda: "../envs/microbeannotator.yml"
    log: "logs/{project}/annotation/microbeAnnotator/{subsample}.log"
    shell:
        """
        mkdir -p {params.dir}
        microbeannotator -i {input.ref} -d {input.db} -o {params.dir} -m {params.method} -t {threads} 2> {log}
        """

rule anvio_script_reformat:
    input: get_protein_fasta
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

rule anvio_run_kegg_kofams:
    input:
        done=rules.anvio_gen_contigs_db.output.done,
        kofam=rules.anvio_setup_kegg_kofams.output
    output:"/tmp/{genome}/anvio_run_kegg_kofams.0"
    params:
       db=rules.anvio_gen_contigs_db.output.db
    conda:"../envs/anvio.yml"
    log: "logs/annotation/anvio_run_kegg_kofams/{genome}.log"
    threads: 40
    shell:
        """
        anvi-run-kegg-kofams -c {params.db} --kegg-data-dir {input.kofam} -T {threads} --just-do-it 2> {log}
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
