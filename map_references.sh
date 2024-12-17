#!/bin/bash

REF_DIR=$1
FASTQ_DIR=$2
OUTPUT_DIR=$3

NEXTFLOW_DIRECTORY=/mnt/extra_disk_03_5TB/JOAO/TOOLS/mapping_workflow


nextflow run $NEXTFLOW_DIRECTORY"/map_references_bash.nf" --references $REF_DIR --reads $FASTQ_DIR --output_dir $OUTPUT_DIR

python $NEXTFLOW_DIRECTORY"/compress_output.py" --output_dir $OUTPUT_DIR


