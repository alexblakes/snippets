# GEL snippets

## Directories

```bash
DIR_AGGV2="/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/genomic_data"
FILE_AGGV2_CHUNK="gel_mainProgramme_aggV2_chr19_648365_2013985.vcf.gz" # Example AggV2 chunk
FILE_GRCH38_FASTA="/public_data_resources/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
FILE_GRCH38_FAI="/public_data_resources/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai"
```

## Job submission
Default LSF job script header. For use with ```bsub < job_script.sh```.

```bash
#!/usr/bin/env bash
#BSUB -q re_gecip
#BSUB -P re_gecip_enhanced_interpretation  
#BSUB -o <path_to/job.%J.%I.out>  
#BSUB -e <path_to/job.%J.%I.err>  
#BSUB -J <jobName>
#BSUB -R "rusage[mem=1G] span[hosts=1]"  
#BSUB -M 1G
#BSUB -n <number_of_cores>  
#BSUB -cwd <"your_dir">

set -euo pipefail

# Set environment
source activate <conda_env> # Not required if submitting from inter queue and conda env already active 
module load <moduleName>

export TMPDIR="/re_scratch/re_gecip/enhanced_interpretation/ablakes/$LSB_JOBNAME"

# My commands...
```

Default LSF bsub command.
Useful for submitting scripts which need command line args.

```bash
bsub \
  -q re_gecip \
  -P re_gecip_enhanced_interpretation \
  -o <path_to/job.%J.%I.out> \
  -e <path_to/job.%J.%I.err> \
  -J <jobName> \
  -R "rusage[mem=1G] span[hosts=1]" \
  -M 1G \
  bash <my_script.sh> <arg1> <arg2> <...>
```

## Using arguments in text files for LSF jobs
Submit a job array, with jobs equal to the number of lines in a file.
(This is in the job submission script.)

```bash
NJOBS=$(wc -l $my_file)

bsub ... -J my_job[1-"${NJOBS}"] ...
```

Get variables from a text file, based on the job index.
(This is in the job script itself.)

```bash
my_var=$(awk -v ix="${LSB_JOBINDEX}" 'NR==ix {print$0}')
```