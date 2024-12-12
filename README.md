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
    v9["minimap2_params"]
    v12["Channel.fromPath"]
    v15["Channel.fromPath"]
    end
    v7([QCReadsProxy])
    v10([MapMinimap2])
    v11([ExtractMappingStatistics])
    v19([CompileMappingStatistics])
    subgraph " "
    v20[" "]
    end
    v1(( ))
    v4(( ))
    v13(( ))
    v18(( ))
    v0 --> v1
    v3 --> v4
    v6 --> v7
    v1 --> v7
    v7 --> v4
    v9 --> v10
    v4 --> v10
    v10 --> v11
    v11 --> v18
    v12 --> v13
    v15 --> v13
    v13 --> v19
    v18 --> v19
    v19 --> v20
```
