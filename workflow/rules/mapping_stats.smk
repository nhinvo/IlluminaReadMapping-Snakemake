rule mapping_stats:
    """
    Obtain mapping stats.
    """
    input: scratch_dict["read_mapping"] / "{sample}_sorted.bam", 
    output: 
        pos = scratch_dict["mapping_stats"] / "{sample}_pos.tsv", 
        contig = scratch_dict["mapping_stats"] / "{sample}_contig.tsv", 
        stats = scratch_dict["mapping_stats"] / "{sample}_mapping_stats.tsv", 
    conda: "../envs/samtools.yaml"
    shell: 
        """
        # obtain depth at all positions 
        samtools depth --threads {resources.cpus_per_task} -o {output.pos} {input}

        # obtain reads mapped to each contig 
        samtools idxstats --threads {resources.cpus_per_task} {input} > {output.contig}

        # obtain mapping stats 
        samtools stats --threads {resources.cpus_per_task} {input} | grep ^SN | cut -f 2- > {output.stats}
        """

rule final:
    """
    Combine all stats. 
    """
    input: 
        pos = expand(scratch_dict["mapping_stats"] / "{sample}_pos.tsv", sample=SAMPLES), 
        contig = expand(scratch_dict["mapping_stats"] / "{sample}_contig.tsv", sample=SAMPLES), 
        stats = expand(scratch_dict["mapping_stats"] / "{sample}_mapping_stats.tsv", sample=SAMPLES), 
    output: 
        final_excel = results_dict["final"],
        summary_mapping_stats = results_dir / "SummaryMapping.png",
        cov_distribution = results_dir / "CovDistribution.png",
    conda: "../envs/data-parse.yaml"
    script: "../scripts/final.py"

