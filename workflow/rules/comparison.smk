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

rule generate_f1_table:
    input:
        config=rules.make_peptable.output,
        force_order=get_genomes,
        metadata=rules.download_gtdb_meta.output
    params:
        genes="results/annotation/anvio/get_sequences_for_gene_calls/",
        genes_kora="results/annotation/kofamscan/functions/default/",
        genes_kore="results/annotation/kofamscan/functions/refined/",
        genes_mara="results/annotation/microbeAnnotator/functions/default/",
        genes_mare="results/annotation/microbeAnnotator/functions/refined/",
        genes_anra="results/annotation/anvio/functions/default/",
        genes_annr="results/annotation/anvio/functions/no_hueristic/",
        genes_anst="results/annotation/anvio/functions/stray/",
        met_kora="results/annotation/kofamscan/metabolism/default/",
        met_mara="results/annotation/microbeAnnotator/metabolism/default/",
        met_mare="results/annotation/microbeAnnotator/metabolism/refined/",
        met_anst="results/annotation/anvio/metabolism/default/",
        met_anra="results/annotation/anvio/metabolism/stray/"
    output: "results/comparison/f1_table/f1_table.tsv"
    log: "logs/comparison/generate_f1_table/log.txt"
    shell:
        """
        sed '1d' config/gtdb_sample_table.tsv > /tmp/gtdb_sample_table.tsv
        echo -e "accession\tgtdb_family\tspecies\tnumber_genes\tgenes_annotated_kora\tgenes_annotated_kore\tgenes_mara\tgenes_mare\tgenes_anra\tgenes_anst\tgenes_annr\tkora_metabolism\tmara_metabolism\tmare_metabolism\tanra_metabolism\tanst_metabolism" > /tmp/f1.table 2> {log} 
        while IFS=$'\t' read -r acc fam _; do
            gene_num=$(grep -c ">" "{params.genes}$acc.faa")
            kofam_genes=$(grep '*' {params.genes_kora}$acc.tsv | cut -f2 | sort -u | wc -l)
            kofam_genes_r=$(grep '*' {params.genes_kore}$acc.tsv | cut -f2 | sort -u | wc -l)
            microbea_genes=$(cat {params.genes_mara}$acc/kofam_results/$acc.faa.kofam.filt | cut -f1 | sort -u | wc -l)
            microbea_genes_r=$(cat {params.genes_mare}$acc/kofam_results/$acc.faa.kofam.filt | cut -f1 | sort -u | wc -l)
   
            species=$(grep $acc {input.metadata} | cut -f2 | sed 's/.*s__//g')
            anra=$(cat {params.genes_anra}$acc.tsv | cut -f1 | sort -u | wc -l)
            anst=$(cat {params.genes_anst}$acc.tsv | cut -f1 | sort -u | wc -l)
            annr=$(cat {params.genes_annr}$acc.tsv | cut -f1 | sort -u | wc -l)
    
            kora_m=$(cut -f10 {params.met_kora}$acc/anvio_estimate_metabolism_modules.txt | awk '{{ sum += $1 }} END {{ print sum }}' | bc -lq)
            mara_m=$(cut -f10 {params.met_mara}$acc/anvio_estimate_metabolism_modules.txt | awk '{{ sum += $1 }} END {{ print sum }}' | bc -lq)
            mare_m=$(cut -f10 {params.met_mare}$acc/anvio_estimate_metabolism_modules.txt | awk '{{ sum += $1 }} END {{ print sum }}' | bc -lq)
            anst_m=$(cut -f10 {params.met_anst}$acc/anvio_estimate_metabolism_modules.txt | awk '{{ sum += $1 }} END {{ print sum }}' | bc -lq)
            anra_m=$(cut -f10 {params.met_anra}$acc/anvio_estimate_metabolism_modules.txt | awk '{{ sum += $1 }} END {{ print sum }}' | bc -lq)
    
            echo -e "$acc\t$fam\t$species\t$gene_num\t$kofam_genes\t$kofam_genes_r\t$microbea_genes\t$microbea_genes_r\t$anra\t$anst\t$annr\t$kora_m\t$mara_m\t$mare_m\t$anra_m\t$anst_m" >> /tmp/f1.table
        done < /tmp/gtdb_sample_table.tsv
        #mkdir -p $(dirname {output})
        grep -v "sample_name" /tmp/f1.table > {output}
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


