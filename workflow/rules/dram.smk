include: "sample_selection.smk"
include: "setup.smk"

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
