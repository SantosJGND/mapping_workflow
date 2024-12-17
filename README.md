# Map References Tool

This tool is designed to map sequencing reads to reference genomes using Nextflow on the server. It leverages conda for environment management and Nextflow for workflow orchestration.

## Prerequisites

Ensure you have the following installed on `genome1`:

- Conda
- Nextflow

## Usage

To use the tool, run the `map_references` script with the following positional arguments:

1. Path to the reference directory
2. Path to the FASTQ reads directory
3. Path to the output directory

### Example

```sh
map_references /path/to/references /path/to/reads /path/to/output
```

#### Pipeline Dependencies

The pipeline requires the following software:

- minimap2
- prinseq++
- samtools

and the pandas python package.

These are included in the conda environment defined in the `environment.yml` file.

#### Pipeline Overview

```mermaid
flowchart TB
    subgraph " "
    v0["Channel.fromPath"]
    v3["Channel.fromPath"]
    v6["prinseq_params"]
    v8["nanofilt_params"]
    v11["minimap2_params"]
    v14["Channel.fromPath"]
    v17["Channel.fromPath"]
    end
    v7([QCReadsPrinseq])
    v9([QCReadsNanofilt])
    v12([MapMinimap2])
    v13([ExtractMappingStatistics])
    v21([CompileMappingStatistics])
    subgraph " "
    v22[" "]
    end
    v1(( ))
    v4(( ))
    v15(( ))
    v20(( ))
    v0 --> v1
    v3 --> v4
    v6 --> v7
    v1 --> v7
    v7 --> v9
    v8 --> v9
    v9 --> v4
    v11 --> v12
    v4 --> v12
    v12 --> v13
    v13 --> v20
    v14 --> v15
    v17 --> v15
    v15 --> v21
    v20 --> v21
    v21 --> v22
```
