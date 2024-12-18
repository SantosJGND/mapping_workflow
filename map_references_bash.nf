#!/usr/bin/env nextflow

/*
* Nextflow script to map a ONT fastq fils to several references after perforing quality control.
*/

/*
*params.output_dir = '/home/bioinf/Desktop/CODE/INSA/TOOLS/map_to_reference/results'
*params.reads = 'data/*.fastq.gz'
*params.references = 'data/references/*.fa'
*/

params.reads = params.reads ?: 'input_reads'
params.references = params.references ?: 'input_references'
params.output_dir = params.output ?: 'output'
params.prinseq_params = '--lc_entropy 0.5 --lc_dust 0.7'
params.minimap2_params = '-ax map-ont --secondary=no'
params.nanofilt_params = '-q 8 -l 50 --headcrop 30 --tailcrop 30 --maxlength 50000'


workflow {

    Channel
        .fromPath("${params.reads}/*.fastq.gz")
        .map { file -> tuple(file.baseName, file) }
        .ifEmpty { error('Cannot find any fastq file') }
        .set { reads_ch }

    // Create a channel for the reference files
    Channel
        .fromPath("${params.references}/*.fa")
        .map { file -> tuple(file.baseName, file) }
        .ifEmpty { error('Cannot find any reference file') }
        .set { reference_ch }


    // QC NANOFILT CHANNEL
    qc_nanofilt_channel = QCReadsNanofilt(reads_ch, params.nanofilt_params)
    // QC PRINSEQ CHANNEL
    qc_channel = QCReadsPrinseq(qc_nanofilt_channel, params.prinseq_params)


    qc_channel = qc_channel.combine(reference_ch)

    mapping_ch = MapMinimap2(qc_channel, params.minimap2_params)
    extract_channel = ExtractMappingStatistics(mapping_ch)

    // COMPILE MAPPING STATISTIC

    CompileMappingStatistics(extract_channel)
}



/*
* Quality control of the reads using nanofilt
*/
process QCReadsNanofilt {
    publishDir "${params.output_dir}/qc_nanofilt_reads", mode: 'copy'

    input:
    tuple val(query_id), path(fastq)
    val nanofilt_params

    output:
    tuple val(query_id), path("${query_id}_good.fastq.gz")

    script:
    """
    zcat ${fastq} | NanoFilt ${nanofilt_params} | bgzip > ${query_id}_good.fastq.gz
    """
}

/*
* Quality control of the reads using prinseq++
*/
process QCReadsPrinseq {
    publishDir "${params.output_dir}/qc_prinseq_reads", mode: 'copy'
    debug true

    input:
    tuple val(query_id), path(fastq)
    val prinseq_params

    output:
    tuple val(query_id), path("${query_id}_good.fastq")

    script:
    """
    prinseq++ ${prinseq_params} -fastq ${fastq} -out_good ${query_id}_good.fastq -out_bad ${query_id}_bad.fastq
    """
}


/*
* QC reads proxy to save time, just copy to same output
*/
process QCReadsProxy {
    tag "QCReadsProxy ${query_id}"

    publishDir "${params.output_dir}/qc_proxy_reads", mode: 'copy'

    input:
    tuple val(query_id), path(fastq)
    val prinseq_params

    output:
    tuple val(query_id), path("${query_id}_good.fastq")

    script:
    """
    cp ${fastq} ${query_id}_good.fastq
    """
}





/*
* Map to a reference using minimap2
*/
process MapMinimap2 {
    tag "MapMinimap2 ${fastq} ${reference_id}"

    publishDir "${params.output_dir}/mapped_reads", mode: 'copy'

    input:
    tuple val(query_id), path(fastq), val(reference_id), path(reference)
    val minimap2_params

    output:
    tuple path("${fastq}_${reference_id}.bam"), val(query_id), val(reference_id)

    script:
    """
    minimap2 ${minimap2_params} ${reference} ${fastq} | samtools view -bS - > ${fastq}_${reference_id}.bam
    """
}
/*
* Map to a reference using bwa
*/
process MapBWAMem {
    tag "MapBWAMem ${fastq} ${reference_id}"

    publishDir "${params.output_dir}/mapped_reads", mode: 'copy'

    input:
    path fastq
    tuple val(reference_id), path(reference)

    output:
    path "${fastq}_${reference_id}.bam"

    script:
    """
    bwa mem -t 4 -R '@RG\\tID:foo\\tSM:bar\\tLB:library1' ${reference} ${fastq} | samtools view -bS - > ${fastq}_${reference_id}.bam
    """
}


/*
* Extract mapping statistics from bam file using samtools
*/
process ExtractMappingStatistics {
    publishDir "${params.output_dir}/mapping_stats", mode: 'copy'

    input:
    tuple path(bam), val(query_id), val(reference_id)

    output:
    tuple path("${bam.baseName}.txt"), val(query_id), val(reference_id)

    script:
    """
    samtools flagstat -O tsv ${bam} > ${bam.baseName}.txt
    """
}


/*
* Extract mapping statistics to single file using python
*/
process CompileMappingStatistics {
    tag "CompileMappingStatistics"

    publishDir "${params.output_dir}/mapping_stats", mode: 'copy'

    input:
    tuple path(file), val(sample_id), val(ref_id)

    output:
    path "${sample_id}_${ref_id}_flagstat.tsv"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd 
    import os
    sample = "${sample_id}"
    reference = "${ref_id}"
    file = "${params.output_dir}/mapping_stats/${sample_id}_good.fastq_${ref_id}.txt"
    output = "${sample_id}_${ref_id}_flagstat.tsv"

    output_df = pd.DataFrame()
    if os.path.exists(output):
        output_df = pd.read_csv(output, sep="\t")


    df = pd.read_csv(file, sep="\t", header=None, names=["value", "qual", "metric"])
    df["sample"] = sample
    df["reference"] = reference

    if output_df.empty:
        output_df = df
    else:
        output_df = pd.concat([output_df, df])

    output_df.to_csv(
        output,
        index=False,
        sep="\t",
    )
    """
}
