"""
Rules for sample selection and generation of the input PEP files needed to run the workflow.
Rules that are for minor file processing steps such as compression are placed directly under
the rule they have input from.

Author: Kathryn Kananen
"""


# Download the gtdb data used for classifying family for each selected sample. From this 
# file we pull the samples being used. 
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

# We supply a file called families_to_sample.txt to randomly select only families we want
# based off of gtdb classification.
rule parse_gtdb_metadata:
    input:
        families=config["families"],
        metadata=rules.untar_metadata.output
    output:temporary("results/temporary/parse_gtdb_metadata/bac120_metadata_r214.families")
    conda:"../envs/sample_selection.yml"
    shell:
        """
        for family in $(cat {input.families})
        do
            grep $family {input.metadata} >> {output}
        done
        """

# Randomly select samples with the use of the shuf command. This selects the minimum amount
# of samples based on the given families presence in the gtdb file. with our chosen
# families we generate a larger sample collection and a smaller one based on 36
rule calculate_subsamples:
    input:
        families=config["families"],
        metadata=rules.parse_gtdb_metadata.output
    output:
        large=temporary("results/temporary/calculate_subsamples/subsamples.large.tsv"),
        set_36=temporary("results/temporary/calculate_subsamples/subsamples.36.tsv")
    conda:"../envs/sample_selection.yml"
    shell:
        """
        n=$(grep -oFf {input.families} {input.metadata} | sort | uniq -c | sort -n | awk 'NR==1 {{print $1}}')
        while read line
        do
            grep $line {input.metadata} > /tmp/calculate_subsamples.$line
            shuf -n $n /tmp/calculate_subsamples.$line >> {output.large}
            shuf -n 36 /tmp/calculate_subsamples.$line >> {output.set_36}
        done < {input.patterns}
        """

# For here we only use the 36 for minimal population size. If a larger analysis is requested
# Then the input here can be linked to the larger sample file. Please note, there may be
# no overlap in the 36 samples selected compared to the larger sample set since the
# samples are randomly selected from the gtdb file independantly. We also keep the 
# genomes that were selected here as a checkpoint to go back to. Users can use the 'touch'
# flag until this point to start the workflow from here.
rule make_genomes_list:
    input: 
        samples=rules.calculate_subsamples.output.set_36,
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

# The pepfile has two columns, the species it is derived from and the name of the genome
# This file needs to be present originally but will be overwritten here and filled with
# the correct subsamples from subsequent steps.
rule make_peptable:
    input:rules.calculate_subsamples.output
    output:"config/gtdb_sample_table.tsv"
    shell:
        """
        printf "%s\t%s\n" "sample_name" "gtdb_family" > {output}
        cut -f17 {input} | sed 's/.*f__//g' | sed 's/;.*//g' > /tmp/make_peptable.family
        cut -f1 {input} | sed 's/.*_G/G/g'> /tmp/make_peptable.genome
        paste /tmp/make_peptable.genome /tmp/make_peptable.family > /tmp/make_peptable.all
        cat /tmp/make_peptable.all >> {output}
        """

# Now we reprocess the DAG to run for all of the subsamples run.
checkpoint process_check:
    input:rules.make_peptable.output
    output:"/tmp/process_check.done"
    shell:"touch {output}"

# Download the genomes from ncbi
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
rule gunzip:
    input:rules.download_genomes.output
    output:"resources/genomes/raw/{genome}.frn"
    conda:"../envs/sample_selection.yml"
    shell:
        """
        gunzip -c {input} > {output}
        """
