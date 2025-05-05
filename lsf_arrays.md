# Idiomatic job arrays on the GEL HPC (LSF)

## Background
I have previously found working with job arrays slightly tedious and intimidating.
I felt that the setup was too involved, and the benefit wasn't worth it.
As a work-around, I had been submitting a few larger jobs to the medium queue, which would complete within a day or so.
However, doing more work with individual-level data across the cohort, I've put together a workflow involving job arrays which is quicker, neater and more reproducible.
I hope it's useful to other GEL users.

## Approach
- Submit an array of, say, 500 jobs (parallelisation #1).
- Within each job process, say, 50 files (parallelisation #2)

## Key themes
- Basic unix commands e.g. `find`, `split`
- Passing arguments to bash scripts and functions
- Bash functions
- GNU parallel
- AWK and the `NR` variable
- LSF variables: `$LSB_JOBINDEX`
- LSF job arrays

## Submission script
We write a script which is used to submit an array job to the GEL HPC.
In this example, we will do some simple filtering on a large number of VCFs. 
We might be doing some analysis on the whole GMS cohort, for example.

The first lines are boilerplate code, including a shebang line and bash "safe mode"
```bash
#!/usr/bin/env bash
set -euo pipefail
```

Next we set some constants which are used throughout the script.
```bash
DIR_VCFS="/path/to/vcf/parent/directory"
DIR_TMP="/path/to/tmp/dir" # Usually in scratch. THIS DIR GETS WIPED. DO NOT USE $pwd
FILE_SPLIT_FILE_PATHS="${DIR_TMP}/split_file_paths.txt" # For example
N=50 # Number of files processed per job
```

To keep things tidy, intermediate outputs will go to the temporary directory described above.
For more security, you can use `mktemp`, but this can be a nuisance to clean up.
Above, we have chosen the name for the temporary directory.
Now, we empty it and start from a clean slate.
```bash
# Save interim files to a temporary directory in scratch
[[ -d $DIR_TMP ]] && rm -rf $DIR_TMP # Safer rm
mkdir -p $DIR_TMP
```

The first step is to get the paths of the target files.
We do this using `find`.
We assume that the file name of each VCF contains the string 'pattern', and has the extension '.vcf.gz'
The paths are split into multiple files, each of which contains N paths.
The split files are given the prefix 'split_vcf_paths'
```bash
find $DIR_VCFS -name *pattern*.vcf.gz -type f \
| split -l $N - "${DIR_TMP}/split_vcf_paths"
```

Now we save the paths of all the split files in another text file.
This text file will be given as an argument to the job script.
```bash
find $DIR_TMP -name split_vcf_paths* -type f > $FILE_SPLIT_FILE_PATHS
```

Finally, we submit an array job to the HPC.
The number of jobs in the array is determined by the number of lines in `$FILE_SPLIT_FILE_PATHS`.
The parameters for the bsub command are just for example.
Submitting this job from the command line (rather than using a script header) allows us to insert variables (such as `njobs`) into the command.

The %100 within the -J flag indicates the number of jobs you want to run concurrently.
Use a max of ~%100 to be fair to other HPC users.

```bash
# Some useful variables for the bsub command
job_name="my_job"
now=$(date -Iseconds)
njobs=$(cat $FILE_SPLIT_FILE_PATHS | wc -l) # Not a useless use of cat

bsub \
    -P re_gecip_enhanced_interpretation \
    -q short \
    -R "rusage[mem=1G]" \
    -M 1G \
    -cwd $(pwd) \
    -o "path/to/logs/${my_job}_${now}_%J_%I.o" \
    -e "path/to/logs/${my_job}_${now}_%J_%I.e" \
    -J "${job_name}[1-${njobs}]%100" \
    bash job_script.sh $FILE_SPLIT_FILE_PATHS
```

And here is the complete submission script:

```bash
#!/usr/bin/env bash
# Array job submission script.
# This script may have rough edges; I haven't tested it in the RE.

set -euo pipefail # Bash "safe mode"

DIR_VCFS="/path/to/vcf/parent/directory"
DIR_TMP="/path/to/tmp/dir" # Usually in scratch. THIS DIR GETS WIPED. DO NOT USE $pwd
FILE_SPLIT_FILE_PATHS="${DIR_TMP}/split_file_paths.txt" # For example
N=50 # Number of files processed per job

# Save interim files to a temporary directory in scratch
[[ -d $DIR_TMP ]] && rm -rf $DIR_TMP # Safer rm
mkdir -p $DIR_TMP

# Find the VCFs we are processing
find $DIR_VCFS -name *pattern*.vcf.gz -type f \
| split -l $N - "${DIR_TMP}/split_vcf_paths"

# Save the paths of the split files to a new file
find $DIR_TMP -name split_vcf_paths* -type f > $FILE_SPLIT_FILE_PATHS

# Some useful variables for the bsub command
job_name="my_job"
now=$(date -Iseconds)
njobs=$(cat $FILE_SPLIT_FILE_PATHS | wc -l) # Not a useless use of cat

# Submit the array job
bsub \
    -P re_gecip_enhanced_interpretation \
    -q short \
    -R "rusage[mem=1G]" \
    -M 1G \
    -cwd $(pwd) \
    -o "path/to/logs/${job_name}_${now}_%J_%I.o" \
    -e "path/to/logs/${job_name}_${now}_%J_%I.e" \
    -J "${job_name}[1-${njobs}]%100" \
    bash job_script.sh $FILE_SPLIT_FILE_PATHS
```

## Job script

This script actually runs the filtering command.
It uses GNU parallel to make use of all of the cores on this worker.

Again, we start with some boilerplate code.
```bash
#!/usr/bin/env bash
set -euo pipefail
```

The file containing the paths of the split files is needed as a command line argument to this script.
```bash
split_file_paths=$1 # set -u will catch if $1 is missing
```

Each job processes the vcf paths in one split file.
This one-liner identifies the split file to be used for this job.
It uses the LSF variable `$LSB_JOBINDEX` and the AWK variable `NR`.
`target_file` is a split file containing N VCF paths.
```bash
target_file=$(< $split_file_paths awk -v ix="$LSB_JOBINDEX" 'NR==ix {print $0}')
```

We write the actual filtering step as a bash function.
For example, here we use bcftools to filter for variants where FILTER=="PASS".
It takes one argument - the path of the VCF to filter.

The function also constructs the output file path.
This could be refactored elsewhere, but it's not a big deal.
```bash
function filterWithBcftools(){
    local file_in=$1 # The VCF to filter, declared as a local variable

    # Get the path to the output file
    local dir_out=$(dirname $file_in)
    local bname=$(basename $file_in .vcf.gz)
    local file_out="${dir_out}/${bname}_filtered.vcf.gz" # For example

    # Run bcftools, for example
    < $file_in bcftools view -i 'FILTER=="PASS"' -o $file_out -W=tbi
}
```

GNU parallel runs in a subshell.
We therefore need to export the function so it is available to the subshells.
```bash
export -f filterWithBcftools
```

Finally, we use parallel to run the filtering step on many cores simultaneously.
We have arranged things so that one VCF is filtered per core.
```bash
parallel --arg-file $target_file --jobs 80% filterWithBcftools {}
```

GNU parallel might not be installed in the GERE by default.
I use it through Conda, for example.

Here is the complete job script.
```bash
#!/usr/bin/env bash
# Job script which runs a bcftools command with GNU parallel.
# bsub is intentionally run on the command line, rather than using a #BSUB header.
# Again, this script hasn't been tested - it may be buggy.

set -euo pipefail # Bash "safe mode"

function filterWithBcftools(){
    local file_in=$1 # The VCF to filter, declared as a local variable

    # Get the path to the output file
    local dir_out=$(dirname $file_in)
    local bname=$(basename $file_in .vcf.gz)
    local file_out="${dir_out}/${bname}_filtered.vcf.gz" # For example

    # Run bcftools, for example
    < $file_in bcftools view -i 'FILTER=="PASS"' -o $file_out -W=tbi
}

export -f filterWithBcftools

split_file_paths=$1 # set -u will catch if $1 is missing


# One-liner to find the specific "split file" processed by this job.
# Recall that each "split file" contains 50 VCF paths - one per line.
target_file=$(< $split_file_paths awk -v ix="$LSB_JOBINDEX" 'NR==ix {print $0}')

parallel --arg-file $target_file --jobs 80% filterWithBcftools {}
```
