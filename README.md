# Human Germline WGS Bioinformatics Workflow
**Andy Ding (andy.ding@sickkids.ca), last updated: March 31, 2022**

Adapted from the GATK Best Practice pipeline (GATK-DRAGEN maximum quality mode) and is based on GRCh38/hg38. To begin this workflow, create a new directory (hereafter referred to as `./`) and enter this directory.

### Preparation
1. Set up Google Cloud [gsutil](https://cloud.google.com/storage/docs/gsutil_install#expandable-1) tool and download GRCh38 reference files

   ```
   wget https://storage.googleapis.com/pub/gsutil.tar.gz
   tar -xvzf gsutil.tar.gz
   mv ~/.boto ~/.boto_old
   echo "[GSUtil]" > ~/.boto
   echo "check_hashes = never" >> ~/.boto
   rm gsutil.tar.gz
   ```
   ```
   mkdir ./hg38_reference
   ./gsutil/gsutil cp gs://gcp-public-data--broad-references/hg38/v0/contamination-resources/1000g/* ./hg38_reference/
   ./gsutil/gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38* ./hg38_reference/
   ./gsutil/gsutil cp gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list ./hg38_reference/
   ./gsutil/gsutil cp gs://gcp-public-data--broad-references/hg38/v0/wgs_coverage_regions.hg38.interval_list ./hg38_reference/
   ./gsutil/gsutil cp gs://gcp-public-data--broad-references/hg38/v0/wgs_evaluation_regions.hg38.interval_list ./hg38_reference/
   ./gsutil/gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz* ./hg38_reference/
   ./gsutil/gsutil cp gs://gcp-public-data--broad-references/hg38/v0/dragen_reference/* ./hg38_reference/
   ./gsutil/gsutil cp gs://gcp-public-data--broad-references/hg38/v0/hg38.* ./hg38_reference/
   ./gsutil/gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.* ./hg38_reference/
   ./gsutil/gsutil cp gs://gcp-public-data--broad-references/hg38/v0/hapmap* ./hg38_reference/
   ./gsutil/gsutil cp gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5* ./hg38_reference/
   ./gsutil/gsutil cp gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1* ./hg38_reference/
   ```
   
2. Install the following packages:

   - [Picard](https://github.com/broadinstitute/picard/releases/latest) (the `jar` file required)

   - [samtools](https://www.htslib.org/download/)

   - [DRAGMAP](https://github.com/Illumina/DRAGMAP) (Use `HAS_GTEST=0 make` instead of `make`)
   
   - [GATK](https://github.com/broadinstitute/gatk/releases)
   
   - [CNVnator](https://github.com/abyzovlab/CNVnator) (Note that the [ROOT package](https://root.cern.ch/) is required. An older version of `samtools` is also required: v1.10 works but v1.14 does not work.)
   
   - `bgzip` and `tabix`. These two tools are typically available by default. If not, they are part of [htslib](https://www.htslib.org/download/).
   
   - [R](https://www.r-project.org/) must be available in the command line.
   
3. Download `resources/hg38_scattered_intervals_40.tar.gz` from this GitHub repository to `./`, then unzip: `tar -xvzf hg38_scattered_intervals_40.tar.gz`.

4. Download `wgs_pipeline_config.txt` from this GitHub repository to `./`. Then specify the **absolute** path to each tool and the folder containing reference files (i.e., `./hg38_reference`). Note that this file is **tab**-separated.


---

### Step 1: FASTQ to gVCF

**Input: FASTQ files of a given sample. Output: gVCF file of this sample. Run this step for each sample.**

**Step-by-step instructions (assume sample ID is `sample001`):**

1. Under `fastq_to_gvcf/` of this GitHub repository, download `single_sample_wgs_pipeline1.pl`, `single_sample_wgs_pipeline2.pl` and `single_sample_wgs_pipeline3.pl` scripts to `./`.

   *Advanced*: The submission scripts (`sub_single_sample_wgs*.pl`) can be found in the `sub/` folder.

2. Place the FASTQ files in the `./fastq/[sample ID]` folder (create this folder if it does not already exist).

   - Name of the FASTQ files **must** follow this format: `[sample ID]_readgroup[number]_R[1/2].fastq.gz`. For example, `sample001_readgroup1_R1.fastq.gz`. 
   
     *Advanced*: If the file names of your FASTQ files do not initially conform to the above format, you may use the `ln -s [old name] [new name]` command to create a shortcut, without renaming your FASTQ files. In this command, `[old name]` is the original file name, and  `[new name]` is the name conforming to the format specified above. For example, `ln -s sample001_S8_L008_R1_001.fastq.gz sample001_readgroup1_R1.fastq.gz`.

   - For paired-end reads, R1 and R2 FASTQ files must both exist. For single-end reads, only the R1 FASTQ file must exist.
   
3. Run `single_sample_wgs_pipeline1.pl`.

   - Edit line 7 of the script (`$path`) to specify the **absolute** path of `./`.
   
   - Create a folder to store temporary files, for example, `./temp_files`. Then edit line 8 of the script (`$temp_path`) to specify the **absolute** path to `./temp_files`.
   
       *Advanced*: If your job scheduler support local scratch space, specify `$temp_path` (line 8 of the script) to this scratch space. This reduces network I/O and may speed up this step.
   
       *Advanced*: On [SickKids Research HPC](https://hpc.ccm.sickkids.ca/), specify `$temp_path` as `/localhd/$ENV{'PBS_JOBID'}`.
   
   - Run the job using `perl single_sample_wgs_pipeline1.pl -s [sample ID]` locally or submit it to HPC cluster using a job scheduler (tested with PBS).
   
       *Advanced*: For PBS, you may use `sub_single_sample_wgs1.pl`, edit line 9 (`$i`) to specify your sample IDs, edit line 16 (`$dir`) to the **absolute** path of `./`. All other parameters (memory, CPU cores, wallclock time and scratch space) are pre-specified. Submit using `perl sub_single_sample_wgs1.pl`.
   
   - For a high-coverage human WGS sample, this step was tested using 100 GB of memory, 64 CPU cores and 500 GB of local scratch space. It takes approximately 18-30 hours to finish.
   
4. Run `single_sample_wgs_pipeline2.pl`.

   - Edit line 8 of the script (`$path`) to specify the **absolute** path of `./`.
   
   - Run the job using `perl single_sample_wgs_pipeline2.pl -s [sample ID] -c [chunk]`. `[chunk]` = 1 to 40 (i.e., run this script 40 times per sample, each time specifying different `[chunk]`). This allows parallelization to speed up the GATK HaplotypeCaller step. 
   
       *Advanced*: For PBS, you may use `sub_single_sample_wgs2.pl`, edit line 9 (`$i`) to specify your sample IDs, edit line 16 (`$dir`) to the **absolute** path of `./`. All other parameters (chunks, memory, CPU cores, wallclock time and scratch space) are pre-specified. Submit using `perl sub_single_sample_wgs2.pl`.
   
   - For a high-coverage human WGS sample, this step was tested using 20 GB of memory and 1 CPU core. It takes up to 4 hours to finish (per chunk).

5. Run `single_sample_wgs_pipeline3.pl`.

   - Edit line 7 of the script (`$path`) to specify the **absolute** path of `./`.
   
   - Run the job using `perl single_sample_wgs_pipeline1.pl -s [sample ID]` locally or submit it to HPC cluster using a job scheduler (tested with PBS).
   
       *Advanced*: For PBS, you may use `sub_single_sample_wgs3.pl`, edit line 9 (`$i`) to specify your sample IDs, edit line 16 (`$dir`) to the **absolute** path of `./`. All other parameters (memory, CPU cores, wallclock time and scratch space) are pre-specified. Submit using `perl sub_single_sample_wgs3.pl`.
       
   - For a high-coverage human WGS sample, this step was tested using 16 GB of memory and 1 CPU core. It takes 1.5-2 hours to finish.
   
6. Examine the outputs:

   - Alignment files (BAM, CRAM) and gVCF file of the given sample can be found at `./single_sample_wgs_final/[sample_ID]`.
   
   - QC files can be found at `./intermediate/[sample_ID]`.
   
---





