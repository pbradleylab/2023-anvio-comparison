""" Samples are subsampled from GTDB by taking the upper bounds of the number
of contigs present in each genome. MAGs will generally have a high number of 
contigs while isolates will have a lower number. To select a wide range of familes
to include in the samples, genomes that had 10 or less contigs were grouped and 
summerized by family and species. Families with the top 20 most total genomes 
where at least 10% were ≤10 contigs and families with the top 20 most total genomes
where under 5% were ≤10 contigs were then taken. 

Author: Patrick J. H. Bradley, Kathryn Kanenen
"""

rule parse_gtdb_metadata:
    input:config["gtdb_metadata"]
    output:
        fc="results/parse_gtdb_metadata/all_family_counts.csv",
        go="results/parse_gtdb_metadata/t20_most_genomes_other.csv"
        glc="results/parse_gtdb_metadata/t20_most_genomes_low_completeness.csv"
    conda:"../envs/sample_selection.yml"
    shell:
        """
        Rscipt parse_gtdb_matadata.R -i {input} -f {output.fc} -l {output.glc} -o {output.go}
        """


