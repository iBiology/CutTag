# CUT&Tag

This pipeline implements the workflow described in 
[CUT&Tag Data Processing and Analysis](https://yezhengstat.github.io/CUTTag_tutorial/index.html). 

## Data pre-processing
Users need to run [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
or other QC tools to make sure the QC meets the minimum requirements described in 
the original paper. If needed, technical replicates or lanes should be merged first.

## Alignment
Reads will be aligned to reference genome and E. coli genome. The former needs 
to be provided using flag `-R` or `--reference` and the corresponding basename 
of bowtie2 index, while the latter needs to be provided via flag `-r` or `--spike` 
along with corresponding basename of bowtie2 index. The following bowtie2 index 
and corresponding species are currently available on the cluster:

| species |                               basename                               |
|---------|:--------------------------------------------------------------------:|
| hg19    |            /project/matzuk/genome/hg19/bowtie2.index/hg19            |
| mm10    |            /project/matzuk/genome/mm10/bowtie2.index/mm10            |
| E. coli |          /project/matzuk/genome/E.coli/bowtie2.index/E.coli          |

In this version, it supports up to 3 different read types: 1 control read (input), 
1 or 2 reads from different conditions. For each read type, both single-end (SE) 
and paired-end (PE) reads are fully supported. The control read can be provided 
using `-i` or `--input` flag along with 1 (SE) or 2 (PE) FASTQ files, the reads 
from first condition can be provided using `-p` or `--ip1` flag along with 1 (SE) 
or 2 (PE) FASTQ files, the reads from second condition can be provided using 
`-q` or `--ip2` flag along with 1 (SE) or 2 (PE) FASTQ files. By default, the 
pipeline will use the basename of r1 for each read type without `.fastq.gz` 
extension as a unique identifier (UID) when there are no UIDs were provided using 
`-I`, `-P`, and `-Q`.

By provide reference genome, E. coli genome, all different types reads, the pipeline 
will automatically align reads to both the reference and E. coli genomes (spike-in).

## Alignment results assessing
Sequencing depth, alignment rate, duplication rate, unique library size, and fragment 
size distribution were assessed and the metrics were saved to CSV files named 
**_alignment.summary.csv_** and _**alignment.fragment.length.csv**_.

## Alignment results filtering
All aligned reads were further filtered by set the minimum quality score to 2, thus 
all alignment results that are below this quality score were eliminated from 
further processing.

## Spike-in calibration
The spike-in calibration is automatically invoked, a scale factor was calculated
using the alignment results from E. coli mapping and a genomic coverage BedGraph
file (_**.bg**_) using this scale factor was also generated.

## Peak calling
Peaks were called by using [SEACR](https://github.com/FredHutch/SEACR/) on 
BedGraph file has normalized fragment counts with E. coli read cont. The
top 1% of peaks (ranked by total signal in each block) were also reported. 
These results can be found in _**.peaks.stringent.bed**_ and 
_**.top.0.01.peaks.stringent.bed**_ files.

## Run the pipeline
The pipeline needs to be wrapper into a PBS job submission script and submit to 
the cluster for processing:
```shell
#!/usr/bin/env bash

#PBS -l nodes=1:ppn=8
#PBS -l vmem=32gb
#PBS -l walltime=02:00:00
#PBS -N CutTag
#PBS -j oe

cd "$PBS_O_WORKDIR" || echo "Failed to enter PBS output work directory!"

source /project/matzuk/software/CutTag/venv/environment.sh

cut-tag \
  -i /project/matzuk/software/CutTag/tests/IgG_rep1_1.fastq.gz \
     /project/matzuk/software/CutTag/tests/IgG_rep1_2.fastq.gz \
  -I IgG \
  -p /project/matzuk/software/CutTag/tests/K27m3_rep1_1.fastq.gz \
     /project/matzuk/software/CutTag/tests/K27m3_rep1_2.fastq.gz \
  -P K27m3 \
  -q /project/matzuk/software/CutTag/tests/K4m3_rep1_1.fastq.gz \
     /project/matzuk/software/CutTag/tests/K4m3_rep1_2.fastq.gz \
  -Q K4m3 \
  -s hg19 \
  -R /project/matzuk/genome/hg19/bowtie2.index/hg19 \
  -r /project/matzuk/genome/E.coli/bowtie2.index/E.coli \
  -o /project/matzuk/software/CutTag/tests/results \
  -n 8
```

The above code snippet can be used as a template to wrapper and 
kickoff the pipeline.

If you want more controls of running the pipeline, you can check 
out its usage:
```shell
$ source /project/matzuk/software/CutTag/venv/environment.sh
$ cut-tag -h
usage: cut-tag [-h] [-i INPUT [INPUT ...]] [-I INPUT_UID] 
               [-p IP1 [IP1 ...]] [-P IP1_UID] 
               [-q [IP2 [IP2 ...]]] [-Q IP2_UID]
               [-s SPECIES] [-R REFERENCE] [-r SPIKE] 
               [-o OUTDIR] [-n CORES] [-S QUALITY_SCORE] [-d] [-D]

cut-tag: a pipeline for processing and analyzing CUT&Tag data.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        Path to input read FASTQ file(s).
  -I INPUT_UID, --input_uid INPUT_UID
                        Unique identifier for input read.
  -p IP1 [IP1 ...], --ip1 IP1 [IP1 ...]
                        Path to IP 1 read FASTQ file(s).
  -P IP1_UID, --ip1_uid IP1_UID
                        Unique identifier for IP 1 read.
  -q [IP2 [IP2 ...]], --ip2 [IP2 [IP2 ...]]
                        Path to IP 2 read FASTQ file(s).
  -Q IP2_UID, --ip2_uid IP2_UID
                        Unique identifier for IP 2 read.
  -s SPECIES, --species SPECIES
                        Name of the species, e.g., hg19, mm10, ... .
  -R REFERENCE, --reference REFERENCE
                        Basename of bowtie2 index for reference genome.
  -r SPIKE, --spike SPIKE
                        Basename of bowtie2 index for spike-in genome.
  -o OUTDIR, --outdir OUTDIR
                        Path to output directory.
  -n CORES, --cores CORES
                        Number of CPUs cores, default: 8.
  -S QUALITY_SCORE, --quality_score QUALITY_SCORE
                        Minimum quality score, default: 2.
  -d, --dryrun          Only print out analytical steps with actually 
                        processing the data.
  -D, --debug           Invoke debug mode that print out more messages 
                        and keep intermedia files.
```

You should alwary issue `source /project/matzuk/software/CutTag/venv/environment.sh` 
to have the pipeline's environment activate and properly set before check its usage 
and use it for data processing.
