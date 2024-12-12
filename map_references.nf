#!/usr/bin/env nextflow

/*
* Nextflow script to map a ONT fastq fils to several references after perforing quality control.
*/

params.output_dir = 'results'
params.reads = 'data/*.fastq.gz'
params.references = 'data/references/*.fa'
params.home = '/home/bioinf/Desktop/CODE/INSA/TOOLS/map_to_reference'
params.prinseq_params = '--lc_entropy 0.5 --lc_dust 0.7'
params.minimap2_params = '-ax map-ont --secondary=no'


workflow {

    // replace reads and references with command line arguments if provided

    Channel
        .fromPath(params.reads)
        .map { file -> tuple(file.baseName, file) }
        .ifEmpty { error('Cannot find any fastq file') }
        .set { reads_ch }

    // create a channel with the reference files. Repeat channel for sample in reads_ch
    Channel
        .fromPath(params.references)
        .map { file -> tuple(file.baseName, file) }
        .ifEmpty { error('Cannot find any reference file') }
        .set { reference_ch }


    qc_channel = QCReadsProxy(reads_ch, params.prinseq_params)

    qc_channel = qc_channel.combine(reference_ch)

    mapping_ch = MapMinimap2(qc_channel, params.minimap2_params)
    extract_channel = ExtractMappingStatistics(mapping_ch)

    // COMPILE MAPPING STATISTICS
    Channel
        .fromPath(params.reads)
        .map { file -> tuple(file.baseName, file) }
        .ifEmpty { error('Cannot find any fastq file') }
        .set { new_reads_ch }

    Channel
        .fromPath(params.references)
        .map { file -> tuple(file.baseName, file) }
        .set { new_reference_ch }

    // combine keep only first element of each tuple
    new_reads_ch
        .combine(new_reference_ch)
        .set { combined_ch }

    CompileMappingStatistics(combined_ch, extract_channel.collect())
}


/*
* Extract mapping statistics to single file using python
*/
process CompileMappingStatistics {
    tag "CompileMappingStatistics"

    publishDir "${params.output_dir}/mapping_stats", mode: 'copy'

    input:
    tuple val(sample_id), val(sample_path), val(ref_id), val(ref_path)
    val mapping_stats

    output:
    path "${sample_id}_${ref_id}_flagstat.txt"

    script:
    """
    python3 ${params.home}/utils/extract_mapping_stats.py ${sample_id} ${ref_id} \
        /${params.output_dir}/mapping_stats/${sample_id}_good.fastq_${ref_id}.txt ${sample_id}_${ref_id}_flagstat.txt
    """
}
/*
* Map to a reference using minimap2
*/
process MapMinimap2 {
    tag "MapMinimap2 ${fastq} ${reference_id}"

    publishDir "${params.output_dir}/mapped_reads", mode: 'copy'

    input:
    tuple path(fastq), val(reference_id), path(reference)
    val minimap2_params

    output:
    path "${fastq}_${reference_id}.bam"

    script:
    """
    minimap2 ${minimap2_params} ${reference} ${fastq} | samtools view -bS - > ${fastq}_${reference_id}.bam
    """
}


/*
* QC reads proxy to save time, just copy to same output
*/
process QCReadsProxy {
    tag "QCReadsProxy ${query_id}"

    publishDir "${params.output_dir}/qc_reads", mode: 'copy'

    input:
    tuple val(query_id), path(fastq)
    val prinseq_params

    output:
    path "${query_id}_good.fastq"

    script:
    """
    cp ${fastq} ${query_id}_good.fastq
    """
}

/*
* Quality control of the reads using prinseq++
*/
process QCReads {
    publishDir "${params.output_dir}/qc_reads", mode: 'copy'
    debug true

    input:
    tuple val(query_id), path(fastq)
    val prinseq_params

    output:
    path "${query_id}_good.fastq"

    script:
    """
    prinseq++ ${prinseq_params} -fastq ${fastq} -out_good ${query_id}_good.fastq -out_bad ${query_id}_bad.fastq
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
    path bam

    output:
    path "${bam.baseName}.txt"

    script:
    """
    samtools flagstat -O tsv ${bam} > ${bam.baseName}.txt
    """
}
