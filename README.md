### ONT Reference Mapping deployment

Map multiple ONT reads to multiple reference genomes using minimap2 and extract mapping statistics.

Currently, the nextflow script needs to be edited to specify the input and output directories, and change the `params` section to specify the parameters for minimap2 and prinseq++ (optional).

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
