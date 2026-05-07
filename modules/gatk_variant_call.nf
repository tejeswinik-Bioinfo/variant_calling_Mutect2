#!/usr/bin/env nextflow

process gatk_mutect2 {

    publishDir "${params.variant_call}/${sample_id}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"
    
    input:
    tuple val (metadata), path (tumor_bam), path (tumor_bai)

    output:
    tuple val (metadata), path ("*somatic*"), emit: vcf
    tuple val (metadata), path ("*somatic*.tbi"), emit: vcf_index
    path("*.stats"), emit: stats
    tuple val (metadata), path ("*f1r2*"), emit: f1r2

    script:
    sample_id = metadata.sampleName

    """
    gatk Mutect2 \
    -R ${params.ref} \
    -I ${tumor_bam} \
    --germline-resource ${params.gNOMAD} \
    --panel-of-normals ${params.PON} \
    --f1r2-tar-gz ${sample_id}_f1r2.tar.gz \
    -O ${sample_id}_somatic.vcf.gz

    """
}

process gatk_mutect2_tumor_normal {

    publishDir "${params.variant_call}/${sample_id}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"
    
    input:
    tuple val (metadata), path (tumor_bam), path (normal_bam)

    output:
    tuple val (metadata), path ("*somatic.vcf.gz"), emit: vcf
    tuple val (metadata), path ("*somatic.vcf.gz.tbi"), emit: vcf_index
    path("*.stats"), emit: stats
    path ("*f1r2*"), emit: f1r2

    script:
    tumor_id = metadata.sampleName
    normal_id = normal_bam.basename

    """
    gatk Mutect2 \
    -R ${params.ref} \
    -I ${tumor_bam} \
    -I ${normal_bam} \
    -normal ${normal_id} \
    --germline-resource ${params.gNOMAD} \
    --panel-of-normals ${params.PON} \
    --f1r2-tar-gz ${tumor_id}_f1r2.tar.gz \
    -O ${tumor_id}_somatic.vcf.gz

    """
}

process gatk_getpileupsummaries{
    publishDir "${params.variant_call}/${sample_id}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (metadata), path (bam_file), path (bai_file)
    output:
    tuple val (metadata), path ("*pileup_summary*")

    script:
    sample_id = metadata.sampleName
    """
    export JAVA_OPTS="-Xmx59G"
    gatk --java-options "-Xmx59G" GetPileupSummaries \
    --verbosity DEBUG \
    -I ${bam_file} \
    -R ${params.ref} \
    -L ${params.common_variants} \
    -V ${params.common_variants} \
    -O ${sample_id}_pileup_summary.table \
    2>&1 | tee -a gatk_debug_live.log
    """
}

process gatk_calculatecontamination{
    publishDir "${params.variant_call}/${sample_id}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (metadata), path (pileup_table)
    output:
    tuple val (metadata), path ("*contamination*")

    script:
    sample_id = metadata.sampleName
    """
    gatk CalculateContamination \
    -I ${pileup_table[0]} \
    -O ${sample_id}_contamination.table
    """

}

process gatk_orientationbias{
    publishDir "${params.variant_call}/${sample_id}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (metadata), path (f1r2_file)
    output:
    tuple val (metadata), path ("*orientation_bias*")

    script:
    sample_id = metadata.sampleName
    """
    gatk LearnReadOrientationModel \
    -I ${f1r2_file} \
    -O ${sample_id}_orientation_bias.tar.gz
    """ 
}

process gatk_filtermutectcalls{
    publishDir "${params.variant_call}/${sample_id}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (metadata), path (raw_vcf), path (orientation_bias), path (contamination_table)

    output:
    tuple val (metadata), path ("*filtered_variants.vcf.gz"), path ("*filtered_variants.vcf.gz.tbi")

    script:
    sample_id = metadata.sampleName
    """
    gatk FilterMutectCalls \
    -V ${raw_vcf[0]} \
    -R ${params.ref} \
    --contamination-table ${contamination_table} \
    --ob-priors ${orientation_bias} \
    -O ${sample_id}_filtered_variants.vcf.gz
    """ 
}

process normalization {
    publishDir "${params.variant_call}/normalized_vcf", mode: "copy"
    conda "conda-forge::gsl bioconda::bcftools=1.23.1"

    input:
    tuple val (metadata), path (filtered_vcf), path (filtered_vcf_index)

    output:
    tuple val (metadata), path ("*norm.vcf.gz"), path ("*norm.vcf.gz.tbi")

    script:
    sample_id = metadata.sampleName
    """
    bcftools norm \
    -m -any \
    -f ${params.ref} \
    -O z \
    -o ${sample_id}_norm.vcf.gz \
     ${filtered_vcf}

    bcftools index -t ${sample_id}_norm.vcf.gz
    """
}



process extract_filtered_variants {
    
    publishDir "${params.variant_call}/filtered_variants/${sample_id}", mode: "copy"
    conda "conda-forge::gsl bioconda::bcftools=1.23.1"

    input:
    tuple val (metadata), path (filtered_vcf), path (filtered_vcf_index)
    output:
    tuple val (metadata), path ("${filtered_vcf.simpleName}_extracted.vcf.gz"), path ("${filtered_vcf.simpleName}_extracted.vcf.gz.tbi")

    script:
    sample_id = metadata.sampleName

    """
    bcftools view -f PASS ${filtered_vcf} -O z -o ${filtered_vcf.simpleName}_extracted.vcf.gz
    bcftools index -t ${filtered_vcf.simpleName}_extracted.vcf.gz
    """

}
