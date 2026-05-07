#!/usr/bin/env nextflow

params.input = ""
params.tumorsampleName = ""
params.normalSampleName = ""
params.tumor_R1 = ""
params.tumor_R2 = ""
params.normal_R1 = ""
params.normal_R2 = ""

include {fastp_trim_reads} from "./modules/fastp"
include {bwa_index; bwa_alignment} from "./modules/bwa"
include {gatk_Mark_Duplicates; gatk_base_recalibrator; gatk_applyBQSR; alignment_metrics; insert_size_metrics} from "./modules/gatk_preprocessing"
include {samtools_faidx; gatk_sequenceDictionary} from "./modules/index_files"
include {gatk_mutect2; gatk_mutect2_tumor_normal; gatk_getpileupsummaries; gatk_calculatecontamination; gatk_orientationbias; gatk_filtermutectcalls; normalization; extract_filtered_variants} from "./modules/gatk_variant_call"
include {gatk_select_variants_SNPs; gatk_select_variants_INDELs; gatk_funcotator; DOWNLOAD_SNPEFF_DB; SNPEFF_ANNOTATE } from "./modules/variant_Annotation"

workflow {
    is_csv = params.input
    is_tumor_only = params.tumor_R1?.trim() && params.tumor_R2?.trim() && !params.normal_R1?.trim()
    is_tumor_normal = params.tumor_R1?.trim() && params.tumor_R2?.trim() && params.normal_R1?.trim() && params.normal_R2?.trim()

    // ========== CSV Input Mode ==========
    if (is_csv) {
        println "CSV input mode: ${params.input}"
        // Determine if tumor-normal mode based on presence of 'normal' type
        csv_file = file(params.input)
        has_normal = csv_file.readLines().drop(1).any { line -> 
            cols = line.split(',')
            cols.size() > 4 && cols[4].trim().toLowerCase() == 'normal'
        }
        println "CSV has normal samples: ${has_normal}"
        
        csv_ch = Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            meta = row.subMap(['sampleName', 'pairedEnd'])
            if (row.type) meta.type = row.type
            tuple(meta, file(row.Read_1), file(row.Read_2))
        }
        //csv_ch.view()
       
        preproc_out = preprocessing(csv_ch)
        variant_calling(preproc_out, has_normal)
    }
    // ========== CLI Tumor-Only Mode ==========
    else if (is_tumor_only) {
        println "CLI tumor-only mode: ${params.tumorsampleName}"
        tumor_meta = [sampleName: params.tumorsampleName ?: "tumor_sample", pairedEnd: true]
        tumor_ch = Channel.of(tuple(tumor_meta, file(params.tumor_R1), file(params.tumor_R2)))

        preproc_out = preprocessing(tumor_ch)
        variant_calling(preproc_out,false)
    }
    // ========== CLI Tumor-Normal Mode ==========
    else if (is_tumor_normal) {
        println "CLI tumor-normal mode: tumor=${params.tumorsampleName}, normal=${params.normalSampleName}"
        tumor_meta = [sampleName: params.tumorsampleName ?: "tumor_sample", pairedEnd: true]
        tumor_ch = Channel.of(tuple(tumor_meta, file(params.tumor_R1), file(params.tumor_R2)))

        normal_meta = [sampleName: params.normalSampleName ?: "normal_sample", pairedEnd: true]
        normal_ch = Channel.of(tuple(normal_meta, file(params.normal_R1), file(params.normal_R2)))

        combined_ch = tumor_ch.concat(normal_ch)
        preproc_out = preprocessing(combined_ch)
        variant_calling(preproc_out,true)
    }
    else {
        error "Please provide either:\n" +
              "  1. CSV: --input <samples.csv> (supports multiple samples, tumor-only or tumor-normal with 'type' column)\n" +
              "  2. Tumor-only: --sampleName <name> --tumor_R1 <R1.fastq.gz> --tumor_R2 <R2.fastq.gz>\n" +
              "  3. Tumor-Normal: --sampleName <name> --normalSampleName <name> --tumor_R1 <R1.fastq.gz> --tumor_R2 <R2.fastq.gz> --normal_R1 <R1.fastq.gz> --normal_R2 <R2.fastq.gz>"
    }
}

workflow preprocessing {
    take:
        sample_ch

    main:

        
        trimmed_reads = fastp_trim_reads(sample_ch)
        bwa_index(params.ref)
        samtools_faidx(params.ref)
        gatk_sequenceDictionary(params.ref)

        aligned = bwa_alignment(fastp_trim_reads.out.trimmed)
        dedup = gatk_Mark_Duplicates(aligned)
        recal = gatk_base_recalibrator(gatk_Mark_Duplicates.out)
        //bqsr_ch = gatk_Mark_Duplicates.out.join(gatk_base_recalibrator.out)
        bqsr_out = gatk_applyBQSR(dedup, recal)
        
        alignment_metrics(bqsr_out)
        insert_size_metrics(bqsr_out)

    emit:
        bqsr_out
}

workflow variant_calling {
    take:
        bqsr_out
        is_tumor_normal

    main:

        split_ch = bqsr_out.multiMap { it -> mutect: it; pileup: it }
        mutect_ch = split_ch.mutect
        pileup_ch = split_ch.pileup

        vcf_ch = Channel.empty()
        pileup_res = Channel.empty()
        f1r2_ch = Channel.empty()

        if (is_tumor_normal) {

            println "Running in tumor-normal mode for variant calling"
            // Split tumor and normal channels
            tumor_sample = mutect_ch.filter { meta, bam, bai -> meta.type == 'tumor' || (!meta.type && meta.sampleName.contains("tumor")) }
            normal_sample = mutect_ch.filter { meta, bam, bai -> meta.type == 'normal' || (!meta.type && meta.sampleName.contains("normal")) }

            // Combine tumor and normal into single tuple
            combined = tumor_sample.combine(normal_sample)
                .map { tumor_tuple, normal_tuple ->
                    def tumor_meta = tumor_tuple[0]
                    def tumor_bam = tumor_tuple[1]
                    def normal_bam = normal_tuple[1]
                    tuple(tumor_meta, tumor_bam, normal_bam)
                }

            mutect2_out = gatk_mutect2_tumor_normal(combined)
            pileup_sample = pileup_ch.filter { meta, bam, bai -> meta.type == 'tumor' || (!meta.type && meta.sampleName.contains("tumor")) }
            pileup_res = gatk_getpileupsummaries(pileup_sample)

            vcf_ch = gatk_mutect2_tumor_normal.out.vcf
            f1r2_ch = gatk_mutect2_tumor_normal.out.f1r2
            
        } else {
            println "Running in tumor-only mode for variant calling"
            vcf_ch = gatk_mutect2(mutect_ch)
            pileup_res = gatk_getpileupsummaries(pileup_ch)

            vcf_ch = gatk_mutect2.out.vcf
            f1r2_ch = gatk_mutect2.out.f1r2            
           
        }

        contamination_out = gatk_calculatecontamination(pileup_res)
        orientation_bias_out = gatk_orientationbias(f1r2_ch)

        filtered_ch = vcf_ch.join(orientation_bias_out).join(contamination_out)
        filtered_out = gatk_filtermutectcalls(filtered_ch)
        norm_vcf = normalization(gatk_filtermutectcalls.out)
        norm_vcf.view()
        
        //snp_indel_ch = gatk_select_variants_SNPs.out.mix(gatk_select_variants_INDELs.out)
        //snp_indel_ch.view()
        extracted_filtered_ch = extract_filtered_variants(norm_vcf)
        //extracted_filtered_ch.view()
        
        DOWNLOAD_SNPEFF_DB()
        SNPEFF_ANNOTATE(extracted_filtered_ch)
        SNPEFF_ANNOTATE.out.ann_vcf.view()

        gatk_select_variants_SNPs(SNPEFF_ANNOTATE.out.ann_vcf)
        gatk_select_variants_INDELs(SNPEFF_ANNOTATE.out.ann_vcf)

    emit:
        annotated_variants = SNPEFF_ANNOTATE.out.ann_vcf


}
    





