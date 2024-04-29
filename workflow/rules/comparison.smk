inlcude: "dram.smk"
include: "anvio.smk"
include: "microbeannotator.smk"
include: "kofamscan.smk"


rule eggnog_mapper:
    input:
        ref=rules.anvi_get_sequences_for_gene_calls.output,
        db=rules.eggnog_db.output.outdir,
        dmnd=rules.eggnog_bacteria_db.output.dmnd,
        gff=rules.anvio_make_gff.output
    output:"results/comparison/eggnog/{genome}.emapper.annotations"
    params:
        gff="results/comparison/eggnog/{genome}",
        tmpdir="results/temporary/comparison/eggnog_mapper/{genome}/"
    log:"logs/comparison/eggnog/{genome}.log"
    conda:"../envs/comparison.yml"
    threads:40
    shell:
        """
        mkdir -p {params.tmpdir}
        emapper.py --itype proteins \
             -i {input.ref}  \
             --tax_scope bacteria \
             --dmnd_db {input.dmnd} \
             --data_dir {input.db} \
             -m diamond \
             --decorate_gff yes \
             -o {params.gff} \
             --override \
             --temp_dir {params.tmpdir} \
             --cpu {threads} 2> {log}
        """

rule generate_f1a:
    input:
        eggnog="",
        linker="",
        strays="",
        mara="",
        mare="",
        anra="",
        anre="",
        annr="",
        kora="",
        kore=""
    output:"results/comparison/f1a.png"
    params:
    conda:"../envs/f1a.yml"
    shell:
        """
        Rscript workflows/scripts/f1a.R
        """
#rule generate_f1b:

#rule generate f1c:

#rule generate 2a:

#rule generate 2b:

#rule generate metabolism_comparison:


