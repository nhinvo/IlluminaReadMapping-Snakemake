# IlluminaReadMapping-Snakemake
bowtie2 read mapping and samtools depth coverage Snakemake pipeline for Illumina short reads. 

## Setup
### 1. Install Snakemake and Conda/Mamba  
Install Snakemake and Conda/Mamba following the instructions at this [link](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#:~:text=for%20installing%20Snakemake.-,Installation%20via%20Conda/Mamba,-This%20is%20the). 

### 2.Set up Snakemake Pipeline
#### 1. Edit Experimental Configurations
Edit **config.yaml** file in the `inputs/` Directory:  

Create **samples.tsv** file in the `inputs/` Directory:  

#### 2. Edit Resource Specifications 
Edit **config.yaml** file in the `profile/` Directory:  

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
**temp.tsv:** tab-separated file with the following columns:  


## Workflow
- QC: bbduk read trimming
- Read Mapping: bowtie2
- Coverage: samtools 