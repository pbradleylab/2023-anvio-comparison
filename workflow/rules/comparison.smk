inlcude: "dram.smk"
include: "anvio.smk"
include: "microbeannotator.smk"
include: "kofamscan.smk"


def get_genomes(wildcards):
    out = []
    for sample in pep.sample_table.sample_name.tolist():
        out.append(rules.anvio_estimate_metabolism.output[0].format(genome=sample))
        out.append(rules.anvio_export_functions_no_huerestic.output[0].format(genome=sample))
        out.append(rules.anvio_export_functions_stray.output[0].format(genome=sample))
        out.append(rules.anvio_export_functions.output[0].format(genome=sample))
        out.append(rules.anvio_estimate_metabolism_stray.output[0].format(genome=sample))
        out.append(rules.kofamscan_refined.output[0].format(genome=sample))
        out.append(rules.kofam_estimate_metabolism.output[0].format(genome=sample))
        out.append(rules.microbeannotator_estimate_metabolism.output[0].format(genome=sample))
        out.append(rules.microbeannotator_estimate_metabolism_refined.output[0].format(genome=sample))
        out.append(rules.eggnog_mapper.output[0].format(genome=sample))
    return out


rule eggnog_mapper:
    input:
        ref=rules.anvio_get_sequences_for_gene_calls.output,
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
