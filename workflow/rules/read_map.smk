rule index_genome:
    "Index reference genome. Credit: Konnor von Emster."
    input: config["input"]["reference genome"], 
    output: touch(scratch_dict["genome_index_done"]), 
    conda: "../envs/bowtie2.yaml"
    shell: "bowtie2-build {input} {input}"


rule map_reads:
    """
    Index reference assembly and map reads to obtain coverage for binning. 

    Credit: Konnor von Emster.
    """
    input: 
        trimmed_r1 = scratch_dict["QC"] / "{sample}_1_trimmed.fastq.gz",
        trimmed_r2 = scratch_dict["QC"] / "{sample}_2_trimmed.fastq.gz",
        ref_genome = config["input"]["reference genome"], 
        genome_index = scratch_dict["genome_index_done"], 
    output: 
        temp(scratch_dict["read_mapping"] / "{sample}.sam"),
    conda: 
        "../envs/bowtie2.yaml"
    shell: 
        """
        # map reads 
        bowtie2 \
            --threads {resources.cpus_per_task} \
            -x {input.ref_genome} \
            -1 {input.trimmed_r1} \
            -2 {input.trimmed_r2} \
            -S {output} 
        """