#!/usr/bin/env nextflow

process gatk_Mark_Duplicates{

    publishDir "${params.outdir}/aligned_reads/${sample_id}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val(metadata), path (aligned_sam)

    output:
    tuple val(metadata), path ("*sorted_dedup*.bam"), path ("*sorted_dedup*.bam.bai"), path ("*sorted_dedup*.bam.sbi")

    script:

    sample_id = metadata.sampleName

    """
    gatk MarkDuplicatesSpark \
    -I $aligned_sam \
    -O ${sample_id}_sorted_dedup_reads.bam
    """

}


process gatk_base_recalibrator{

    publishDir "${params.outdir}/aligned_reads/${sample_id}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val(metadata), path (dedup_bam), path (bai), path (sbi)
    output:
    tuple val(metadata), path ("*_recal_data.table")

    script:
    sample_id = metadata.sampleName

    """
    gatk BaseRecalibrator \
    -I ${dedup_bam} \
    -R ${params.ref} \
    --known-sites ${params.known_sites} \
    -O ${sample_id}_recal_data.table
    """

}

process gatk_applyBQSR{
    publishDir "${params.outdir}/aligned_reads/${sample_id}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"
    
    input:
    tuple val(metadata), path (dedup_bam), path (bai), path (sbi)
    tuple val(metadata), path (recal_table)

    output:
    tuple val(metadata), path ("*dedup_bqsr*.bam"), path ("*dedup_bqsr*.bai")

    script:
    sample_id = metadata.sampleName

    """
    gatk ApplyBQSR \
    -I ${dedup_bam} \
    -R ${params.ref} \
    --bqsr-recal-file ${recal_table} \
    -O ${sample_id}_dedup_bqsr.bam
    """ 
}




process alignment_metrics{
    publishDir "${params.outdir}/aligned_reads/${sample_id}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val(metadata), path (dedup_bqsr_bam), path (bai)
    output:
    tuple val(metadata), path ("*alignment_summary*")

    script:
    sample_id = metadata.sampleName

    """
    gatk CollectAlignmentSummaryMetrics \
    -I ${dedup_bqsr_bam[0]} \
    -R ${params.ref} \
    -O ${sample_id}_alignment_summary_metrics.txt
    """

}

process insert_size_metrics{
    publishDir "${params.outdir}/aligned_reads/${sample_id}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val(metadata), path (dedup_bqsr_bam), path (bai)
    output:
    tuple val(metadata), path ("*insert_size_metrics*")

    script:
    sample_id = metadata.sampleName

    """
    gatk CollectInsertSizeMetrics \
    -I ${dedup_bqsr_bam[0]} \
    -O ${sample_id}_insert_size_metrics.txt \
    -H ${sample_id}_insert_size_histogram.pdf
    """


}