#!/bin/bash
#$ -q som,bio
#$ -M freese@uci.edu
#$ -m ea
#$ -cwd
#$ -j y
#$ -N fastq_bam
#$ -pe one-node-mpi 16
#$ -R y

module load picard/2.18.4

f=${p}.fa.gz
o=${p}.bam

set -x

echo $p
echo $f
echo $o

picard FastqToSam \
    F1=$f \
    O=$o \
    SM=$p