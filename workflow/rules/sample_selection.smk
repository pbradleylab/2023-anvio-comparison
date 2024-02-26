""" Samples are subsampled from GTDB by taking the upper bounds of the number
of contigs present in each genome. MAGs will generally have a high number of 
contigs while isolates will have a lower number. To select a wide range of familes
to include in the samples, genomes that had 10 or less contigs were grouped and 
summerized by family and species. Families with the top 20 most total genomes 
where at least 10% were ≤10 contigs and families with the top 20 most total genomes
where under 5% were ≤10 contigs were then taken. 

Author: Patrick J. H. Bradley, Kathryn Kanenen
"""

rule download_metadata:
    output:
        tsv="resources/bac120_metadata_r214.tsv",
        tar="resources/bac120_metadata_r214.tar.gz"
    shell:
        """
        wget -c https://data.gtdb.ecogenomic.org/releases/release214/214.0/bac120_metadata_r214.tar.gz -P $(dirname {output.tsv})
        tar -xzvf {output.tar} -C $(dirname {output.tsv})
        """

rule parse_gtdb_metadata:
    input:rules.download_metadata.output.tsv
    output:
        fc="results/parse_gtdb_metadata/all_family_counts.csv",
        go="results/parse_gtdb_metadata/t20_most_genomes_other.csv",
        glc="results/parse_gtdb_metadata/t20_most_genomes_low_completeness.csv"
    conda:"../envs/sample_selection.yml"
    shell:
        """
        Rscript workflow/scripts/parse_gtdb_metadata.R -i "{input}" -c "{output.fc}" -l "{output.glc}" -o "{output.go}"
        """

rule fetch_genomes_high:
    input: 
        go=rules.parse_gtdb_metadata.output.go,
        meta=rules.download_metadata.output.tsv
    output:"resources/genomes.high"
    shell:
        """
        while read line
        do
            family=$(echo $line | cut -d',' -f5)
            grep $family {input.meta} >> /tmp/fetch_genomes_high.family
        done < {input.go}
        cat /tmp/fetch_genomes.family | cut -f1 | sed 's/.*_G/G/g' >> /tmp/fetch_genomes_high.genomes
        sort -u /tmp/fetch_genomes_high.genomes > {output}
        """

rule fetch_genomes_low:
    input:
        glc=rules.parse_gtdb_metadata.output.glc,
        meta=rules.download_metadata.output.tsv
    output:"resources/genomes.low"
    shell:
        """
        while read line
        do
            family=$(echo $line | cut -d',' -f5)
            grep $family {input.meta} >> /tmp/fetch_genomes_low.family
        done < {input.glc}
        cat /tmp/fetch_genomes.family | cut -f1 | sed 's/.*_G/G/g' >> /tmp/fetch_genomes_low.genomes
        sort -u /tmp/fetch_genomes_low.genomes > {output}
        """

rule download_genomes_refseq:
    input:
        low=rules.fetch_genomes_high.output,
        high=rules.fetch_genomes_low.output
    output:directory("resources/genomes/refseq/")
    conda:"../envs/sample_selection.yml"
    shell:
        """
        cat {input.low} {input.high} > /tmp/genomes.all
        grep "GCF" /tmp/genomes.all > /tmp/genomes.refseq
        ncbi-genome-download -A /tmp/genomes.refseq --section refseq bacteria -F "fasta,protein-fasta" -o $(dirname {output})
        """

rule download_genomes_genbank:
    input:
        low=rules.fetch_genomes_high.output,
        high=rules.fetch_genomes_low.output
    output:directory("resources/genomes/genbank/")
    conda:"../envs/sample_selection.yml"
    shell:
        """
        cat {input.low} {input.high} > /tmp/genomes.all
        grep "GCA" /tmp/genomes.all > /tmp/genomes.genbank
        ncbi-genome-download -A /tmp/genomes.genbank --section genbank bacteria -F "fasta,protein-fasta" -o $(dirname {output})
        """

rule symlink_genomes:
    input:
        genbank=rules.download_genomes_genbank.output,
        refseq=rules.download_genomes_refseq.output
    output:directory("resources/genomes/symlink")
    shell:
        """
        ln -s $PWD/{input.genbank} $PWD/{output}
        ln -s $PWD/{input.refseq} $PWD/{output}
        """

