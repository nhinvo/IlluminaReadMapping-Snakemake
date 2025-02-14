rule samtools_depth:
    """
    Obtain coverage at each position. 
    """
    input: scratch_dict["read_mapping"] / "{sample}_sorted.bam", 
    output: scratch_dict["mapping_stats"] / "{sample}_pos.tsv", 
    conda: "../envs/samtools.yaml"
    shell: "samtools depth --threads {resources.cpus_per_task} -o {output} {input}"


rule samtools_contig_depth:
    """
    Obtain coverage of each contig (useful for mapping to concat genome). 
    """
    input: scratch_dict["read_mapping"] / "{sample}_sorted.bam", 
    output: scratch_dict["mapping_stats"] / "{sample}_contig.tsv", 
    conda: "../envs/samtools.yaml"
    shell: "samtools idxstats --threads {resources.cpus_per_task} {input} > {output}"


rule samtools_stats:
    """
    Obtain summary stats (e.g. total reads, reads mapped, etc.)
    """
    input: scratch_dict["read_mapping"] / "{sample}_sorted.bam", 
    output: scratch_dict["mapping_stats"] / "{sample}_mapping_stats.tsv", 
    conda: "../envs/samtools.yaml"
    shell: "samtools stats --threads {resources.cpus_per_task} {input} | grep ^SN | cut -f 2- > {output}"


rule final:
    """
    Combine all stats. 
    """
    input: 
        pos = expand(scratch_dict["mapping_stats"] / "{sample}_pos.tsv", sample=SAMPLES), 
        contig = expand(scratch_dict["mapping_stats"] / "{sample}_contig.tsv", sample=SAMPLES), 
        stats = expand(scratch_dict["mapping_stats"] / "{sample}_mapping_stats.tsv", sample=SAMPLES), 
    output: results_dict["final"],
    conda: "../envs/data-parse.yaml"
    script: "../scripts/final.py"

