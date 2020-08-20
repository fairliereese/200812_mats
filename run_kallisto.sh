#!/bin/bash
#$ -q som,bio
#$ -M freese@uci.edu
#$ -m ea
#$ -cwd
#$ -j y
#$ -N kallisto
#$ -pe one-node-mpi 16
#$ -R y

module load kallisto/0.46.1

idx=gencode.v29.transcripts.idx
fa=${p}.fa.gz

kallisto quant \
	-i $idx \
	-t 16 \
	-o kallisto_${p} \
	--single \
	-l 100 \
	-s 30 \
	$fa