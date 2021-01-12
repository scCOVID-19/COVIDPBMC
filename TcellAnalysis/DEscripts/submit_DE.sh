#! /usr/bin/bash

## submit the DE testing to the cluster

## Run mylo on a subset of the data
BASE_DIR="/hps/research/sds/sds-marioni-cov/DifferentialExpression"
SRC=$(echo $BASE_DIR"/src")

OUTLOG=$(eval 'echo "$BASE_DIR"/logs/Tcell_DE_testing.out')
ERROR=$(eval 'echo "$BASE_DIR"/logs/Tcell_DE_testing.err')

module load singularity/3.5.0;

cd $SRC

JOB="bsub -q research-rh74 -M 160000 -R "rusage[mem=40000]" -J Tcell_DE_testing -n 1 -T 1500 -o $OUTLOG -e $ERROR singularity exec  -B /hps/research/sds/sds-marioni-cov /hps/research/sds/sds-marioni-cov/singularity/CovidPBMC.sif Rscript $SRC/Compute_CellType_DE_Trend.R"

echo $JOB
eval $JOB

