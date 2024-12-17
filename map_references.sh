#!/bin/bash

REF_DIR=$1
FASTQ_DIR=$2
OUTPUT_DIR=$3

NEXTFLOW_DIRECTORY=/path/to/nextflow

conda activate map_reference

nextflow run $NEXTFLOW_DIRECTORY"/map_references.nf" --references $REF_DIR --reads $FASTQ_DIR --output_dir $OUTPUT_DIR

python $NEXTFLOW_DIRECTORY"compress_output.py" $OUTPUT_DIR