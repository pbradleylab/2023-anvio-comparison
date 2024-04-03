include: "sample_selection.smk"
include: "setup.smk"
import glob

# Because kofamscan is too slow in processing large files, these files need to be split up into smaller chunks ~500Kb to finish in
# in under 3 weeks time. Here we evaluate the sequence and then split it apart to a smaller size. We then recalculate the input and
# add the additional files to the input required of the rule. In post, we then merge the files htat were originally split.

rule split_large_fasta:
    input:rules.transeq.output
    output:directory("results/split_large_fasta/{genome}/")
    conda: "../envs/kofamscan.yml"
    shell:
        """
        size=$(stat -c %s {input})
        if [ $size -gt 512000 ]; then
            seqkit split {input} --by-size 5 -O {output}
        else
            mkdir -p {output}
            ln -s $PWD/{input} {output}/$(basename {input})
        fi
        """

rule kofamscan:
    input:
        ko_list=rules.untar_kofams_list.output,
        profiles=rules.untar_kofams_profile.output,
        indir=rules.split_large_fasta.output
    output:"results/annotation/kofamscan/{genome}.tsv"
    conda:"../envs/kofamscan.yml"
    shell:
        """
        for f in $(ls {input.indir})
        do
            name=$(echo $f | sed 's/.faa//g')
            ls {output} || exec_annotation -f detail-tsv -p {input.profiles} -o /tmp/{wildcards.genome}/$name.out {input.indir}/$f -k {input.ko_list} --tmp-dir /tmp/{wildcards.genome}/$name/
        done
        cat /tmp/{wildcards.genome}/$name.out | uniq > {output}
        """



rule kofamscan_refined:
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
        anvi-estimate-metabolism -c {params.db} --user-modules {input.enzymes} --kegg-data-dir /home/kananen.13/workflows/2023-anvio-comparison/resources/feb_anvio/stray/ -O {params.prefix} 2> {log}
        """
