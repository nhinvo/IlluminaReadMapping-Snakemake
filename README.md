# IlluminaReadMapping-Snakemake
Snakemake Pipeline for Coverage QC of Illumina short reads. 

## Setup
### 1. Install Snakemake and Conda/Mamba  
Install Snakemake and Conda/Mamba following the instructions at this [link](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#:~:text=for%20installing%20Snakemake.-,Installation%20via%20Conda/Mamba,-This%20is%20the). 


### 2.Set up Snakemake Pipeline
#### 1. Edit Experimental Configurations
Edit **config.yaml** file in the `inputs/` Directory:  
  - If needed, edit paths to: "sample table", "reference genome", "scratch directory" and "results directory"

Create **samples.tsv** file in the `inputs/` Directory:  
  - Required columns: 
    - sample: unique sample name 
    - forward read: path to forward (1) raw read file. 
    - reverse read: path to reverse (2) raw read file. 
  - Add as many additional metadata column as needed (these will not affect the pipeline)


#### 2. Edit Resource Specifications 
Edit **config.yaml** file in the `profile/` Directory:  
  - Edit number of jobs to run at once to preferred number at "jobs"
  - Edit "partition" to name of partition you would like to run jobs on 
  - Edit "time", "mem", and "cpus_per_task" to match with your cluster specifications

## Running Pipeline 
### 1. Edit main "run_mapping.sbatch" Script
- Edit name of partition on line 4
- Edit name of conda environment with Snakemake installed on line 10 (if env name is other than "snakemake")

### 2. Submit Script to Cluster
- Submit job to cluster by:  
  ```
  sbatch run_assembly.sbatch
  ```  

## Troubleshooting 
As the pipeline runs, log messages will be saved into file named "main.[slurm_job_ID].err" in the "logs" folder. Here are some tips on debugging: 
- Each rule/job in the pipeline will get its own .err log file. When a rule fails, check the subfolder with the rule's name and sample where it failed. 
- Locked directory: if you get the error message "Error: Directory cannot be locked":
  - Make sure that all jobs from your previous run of the pipeline have completed/cancelled
  - Uncomment line 10 "unlock: True" in `profile/config.yaml` file
  - Run the pipeline following previous step: [Submit Script to Cluster](#2-submit-script-to-cluster)
  - Wait for snakemake to complete running
  - Comment line 10
  - Re-run pipeline, the lock should now be removed  


## Results
**MappingCov.xlsx:** Excel file with 2 sheets of mapping statistics 
1. mapping_stats: summary mapping statistics for each sample. 
2. contig_mapping_stats: percent of reads mapped to each contig for each sample. Useful when mapping to a concatenated genome. 


## Workflow
- QC: bbduk read trimming
- Read Mapping: bowtie2
- Coverage: samtools (depth, idxstats, stats)