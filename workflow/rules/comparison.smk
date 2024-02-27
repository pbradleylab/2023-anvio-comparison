include: "annotation.smk"

rule microbeannotator_clean:
    input:rules.microbeannotator.output
    output:"results/comparissons/microbeannotator/{genome}.clean"
    shell:
        """
        sed '1,3d' {input} | cut -f3 | sort -u > {output}
        """

rule anvio_clean:
    input:rules.anvio_export_functions.output
    output:"results/comparissons/microbeannotator/{genome}.clean"
    shell:
        """
        sed '1d' {input} | cut -f3 | sort -u > {output}
        """

rule microbeA_anvio_diff:
    input:
        microbeannotator=rules.microbeannotator_clean.output,
        anvio=rules.anvio_clean.output
    output:"results/comparissons/microbeA_anvio_diff/{genome}.diff"
    shell:
        """
        grep -F -v -x -f {input.microbeannotator} {input.anvio} > {output}
        """

rule microbeA_kofam_diff:
    input:
        microbeannotator=rules.microbeannotator_clean.output,
        kofam=rules.kofamscan.output
    output:"results/comparissons/microbeA_kofam_diff/{genome}.diff"
    shell:
        """
        grep -F -v -x -f {input.microbeannotator} {input.kofam} > {output}
        """

rule anvio_microbeA_diff:
    input:
        microbeannotator=rules.microbeannotator_clean.output,
        anvio=rules.anvio_clean.output
    output:"results/comparissons/anvio_microbeA_diff/{genome}.diff"
    shell:
        """
        grep -F -v -x -f {input.anvio} {input.microbeannotator} > {output}
        """

rule anvio_kofam_diff:
    input:
        anvio=rules.anvio_clean.output,
        kofam=rules.kofamscan.output
    output:"results/comparissons/anvio_kofam_diff/{genome}.diff"
    shell:
        """
        grep -F -v -x -f {input.anvio} {input.kofam} > {output}
        """

rule kofam_microbeA_diff:
    input:
        kofam=rules.kofamscan.output,
        microbeannotator=rules.microbeannotator_clean.output
    output:"results/comparissons/kofam_microbeA_diff/{genome}.diff"
    shell:
        """
        grep -F -v -x -f {input.kofam} {input.microbeannotator} > {output}
        """

rule kofam_anvio_diff:
    input:
        kofam=rules.kofamscan.output,
        anvio=rules.anvio_clean.output
    output:"results/comparissons/kofam_anvio_diff/{genome}.diff"
    shell:
        """
        grep -F -v -x -f {input.kofam} {input.anvio} > {output}
        """

# Generate a venn diagram of all the different ko reported between 
# the tools by default regradless of their e-value or bitscore.
rule generate_venn_diff_all:
    input:
        rules.kofam_microbeA_diff.output,
        rules.kofam_anvio_diff.output,
        rules.microbeA_anvio_diff.output,
        rules.microbeA_kofam_diff.output,
        rules.anvio_microbeA_diff.output,
        rules.anvio_kofam_diff.output
    output:"results/comparison/generate_venn/kegg_comparison.pdf"
    shell:
        """
        Rscipt venn_generation.R -km -ka -ma -mk -am -ak -o {output}
        """
# Generate a scatterplot of the different bitscores of hits that are 
# found in common accross the groups.
#rule generate_scatter_same_eval:

