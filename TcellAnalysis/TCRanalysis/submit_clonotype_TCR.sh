#! /usr/bin/bash

## Run mylo on a subset of the data

BASE_DIR="/hps/research/sds/sds-marioni-cov"
SRC=$(echo $BASE_DIR"/scripts")

OUTLOG=$(eval 'echo "$BASE_DIR"/logs/TCR_clonotypes.out')
ERROR=$(eval 'echo "$BASE_DIR"/logs/TCR_clonotypes.err')

JOB="bsub -q research-rh74 -M 160000 -R "rusage[mem=20000]" -J TCR_clonotypes -n 1 -T 1500 -o $OUTLOG -e $ERROR Rscript $SRC/tcr_clonotypes.R"

echo $JOB
eval $JOB
