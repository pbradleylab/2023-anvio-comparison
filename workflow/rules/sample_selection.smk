rule download_metadata:
    output:temp("resources/bac120_metadata_r214.tar.gz")
    shell:
        """
        wget -c https://data.gtdb.ecogenomic.org/releases/release214/214.0/bac120_metadata_r214.tar.gz -P $(dirname {output})
        """

rule untar_metadata:
    input:rules.download_metadata.output
    output:"resources/bac120_metadata_r214.tsv"
    shell:
        """
        tar -xzvf {input} -C $(dirname {output})
        """

rule parse_gtdb_metadata:
    input:
        families="config/families_to_sample.txt",
        metadata=rules.untar_metadata.output
    output:"results/temporary/parse_gtdb_metadata/bac120_metadata_r214.families"
    conda:"../envs/sample_selection.yml"
    shell:
        """
        for family in $(cat {input.families})
        do
            grep $family {input.metadata} >> {output}
        done
        """

rule calculate_subsamples:
    input:
        metadata=rules.parse_gtdb_metadata.output,
        patterns="config/families_to_sample.txt"
    output:"results/temporary/calculate_subsamples/subsamples.tsv"
    conda:"../envs/sample_selection.yml"
    shell:
        """
        n=$(grep -oFf {input.patterns} {input.metadata} | sort | uniq -c | sort -n | awk 'NR==1 {{print $1}}')
        while read line
        do
            grep $line {input.metadata} > /tmp/calculate_subsamples.$line
            shuf -n $n /tmp/calculate_subsamples.$line >> {output}
        done < {input.patterns}
        """

rule make_genomes_list:
    input: 
        samples=rules.calculate_subsamples.output,
        meta=rules.untar_metadata.output
    output:"resources/make_genomes_list/genomes.selected"
    shell:
        """
        while read line
        do
            genome=$(echo $line | cut -d' ' -f1 | sed 's/.*_G/G/g')
            echo $genome >> {output}
        done < {input.samples}
        """

rule make_peptable:
    input:rules.calculate_subsamples.output
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
    output:"resources/genomes/raw/{genome}.frn.gz"
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

rule download_gff:
    input:rules.make_genomes_list.output
    output:"resources/genomes/gff/{genome}.gff.gz"
    conda:"../envs/sample_selection.yml"
    shell:
        """
        if [[ "{wildcards.genome}" == "GCA"* ]]; then
            ncbi-genome-download -A "{wildcards.genome}" --section genbank bacteria -F "gff" -o /tmp/{wildcards.genome}
            mv /tmp/{wildcards.genome}/genbank/bacteria/{wildcards.genome}/*.gz {output}
        else
            ncbi-genome-download -A "{wildcards.genome}" --section refseq bacteria -F "gff" -o /tmp/{wildcards.genome}
            mv /tmp/{wildcards.genome}/refseq/bacteria/{wildcards.genome}/*.gz {output}
        fi
        """

rule gunzip:
    input:rules.download_genomes.output
    output:"resources/genomes/raw/{genome}.frn"
    conda:"../envs/sample_selection.yml"
    shell:
        """
        gunzip -c {input} > {output}
        """
