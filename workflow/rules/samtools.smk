rule samtools_view:
    """
    Convert sam to bam.
    """
    input: scratch_dict["read_mapping"] / "{sample}.sam", 
    output: temp(scratch_dict["read_mapping"] / "{sample}_unsorted.bam"), 
    conda: "../envs/samtools.yaml"
    shell:
        """
        samtools view \
            --bam \
            --threads {resources.cpus_per_task} \
            --output {output} \
            {input}
        """


rule samtools_sort_index:
    """
    Sort index bam. 
    """
    input: 
        scratch_dict["read_mapping"] / "{sample}_unsorted.bam", 
    output: 
        sorted_bam = scratch_dict["read_mapping"] / "{sample}_sorted.bam", 
        indexed_bam = scratch_dict["read_mapping"] / "{sample}_sorted.bam.bai", 
    conda: 
        "../envs/samtools.yaml"
    shell:
        """
        # sort bam 
        samtools sort \
            --threads {resources.cpus_per_task} \
            -o {output.sorted_bam} \
            {input}

        # index the sorted_bam
        samtools index \
            --threads {resources.cpus_per_task} \
            --bai \
            --output {output.indexed_bam} \
            {output.sorted_bam}
        """