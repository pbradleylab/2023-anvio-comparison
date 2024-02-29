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

rule calculate_subsamples:
    input:
        low=rules.parse_gtdb_metadata.output.glc,
        high=rules.parse_gtdb_metadata.output.go
    output:"results/genomes/calculate_subsamples/subsamples.tsv"
    conda:"../envs/sample_selection.yml"
    shell:
        """
        cat {input.low} > /tmp/genomes.all
        sed '1d' {input.high} >> /tmp/genomes.all
        python3 workflow/scripts/subsample_gtdb.py /tmp/genomes.all {output}
        """

rule make_genomes_list:
    input: 
        samples=rules.calculate_subsamples.output,
        meta=rules.download_metadata.output.tsv
    output:"resources/genomes.selected"
    shell:
        """
        while read line
        do
            family=$(echo $line | cut -d',' -f1)
            grep $family {input.meta} >> /tmp/$family.make_genomes_list
        done < {input.samples}

        while read line
        do
            n=$(echo $line | cut -d',' -f3)
            family=$(echo $line | cut -d',' -f1)
            shuf -n $n /tmp/$family.make_genomes_list >> {output}
        done < {input.samples}
        """

rule make_peptable:
    input:rules.make_genomes_list.output
    output:"config/gtdb_sample_table.tsv"
    shell:
        """
        echo "sample_name\tfamily" > {output}
        cut -f17 {input} | sed 's/.*f__//g' | sed 's/;.*//g' > /tmp/make_peptable.family
        cut -f1 {input} | sed 's/.*_G/G/g'> /tmp/make_peptable.genome
        paste /tmp/make_peptable.genome /tmp/make_peptable.family > /tmp/make_peptable.all
        cat /tmp/make_peptable.all >> {output}
        """

checkpoint process_check:
    input:rules.make_peptable.output
    output:"/tmp/process_check.done"
    shell:"touch {output}"

rule download_genomes:
    input:rules.make_genomes_list.output
    output:"resources/genomes/{genome}.frn.gz"
    conda:"../envs/sample_selection.yml"
    shell:
        """
        if [[ "{wildcards.genome}" == "GCA"* ]]; then
            ncbi-genome-download -A "{wildcards.genome}" --section genbank bacteria -F "fasta" -o /tmp/{wildcards.genome}
            mv /tmp/{wildcards.genome}/genbank/bacteria/{wildcards.genome}/*.gz {output}
        else
            ncbi-genome-download -A "{wildcards.genome}" --section refseq bacteria -F "fasta" -o /tmp/{wildcards.genome}
            mv /tmp/{wildcards.genome}/refseq/bacteria/{wildcards.genome}/*.gz {output}
        fi
        """

rule gunzip:
    input:rules.download_genomes.output
    output:"resources/genomes/{genome}.frn"
    conda:"../envs/sample_selection.yml"
    shell:
        """
        gunzip -c {input} > {output}
        """

rule transeq:
    input:rules.gunzip.output
    output:"resources/genomes/symlink/{genome}/{genome}.faa",
    conda:"../envs/sample_selection.yml"
    shell:
        """
        transeq {input} {output} 
        """
