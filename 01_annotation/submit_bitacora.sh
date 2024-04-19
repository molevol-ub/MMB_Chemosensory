#!/bin/bash
#$ -cwd
#$ -V
#$ -N 'bitacora'
#$ -q h13.q
#$ -pe ompi255h13 8

source ~/.bash_profile
#conda activate minibusco

DB=GITHUB_DIR/01_annotation/DB
GENOME=GITHUB_DIR/01_annotation/Files/Sinv_genome.softMasked.fasta
GFF=GITHUB_DIR/01_annotation/Files/Sinv_annotation.gff3
PROT=GITHUB_DIR/01_annotation/Files/Sinv_annotation.pep_noiso.fasta

perl GITHUB_DIR/programs/bitacora_1.4.1/runBITACORA_command_line.sh -m full -q $DB -g $GENOME -f $GFF -p $PROT -n Sinv -t 8 -r F

