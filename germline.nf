#!/usr/bin/env nextflow
nextflow.enable.dsl=2

VERSION="0.0"

log.info "===================================================================="
log.info "FASTQC fasta files in folder                                        "
log.info "===================================================================="

params.help = ""
if (params.help) {
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "nextflow run step_A.nf --fasta list.txt"
  log.info " "
    log.info "   list.txt contains the name with full path of a set of samples"
    log.info "   only .gz files are supported"
    log.info " "
  log.info "===================================================================="
  exit 1
}


//includeConfig 'step_A.config'
include {fastqc; multiqc } from './modules/qc.nf'
include {trim_illumina } from './modules/trim.nf'
include {split_align_mark; alignment_metrics; base_recalibrator } from './modules/gatk.nf'
include {haplotype_caller_gvcf}  from './modules/gatk.nf'
include {haplotype_caller_gvcf; genomicsDBImport } from './modules/gatk.nf'
include {gather; variantRecalibratorSNP; variantRecalibratorINDEL } from './modules/gatk.nf'
include {applyVqsrINDEL; applyVqsrSNP } from './modules/gatk.nf'
include {vep_annotation } from './modules/annotation.nf'


workflow annotate {

    take:
    vcf

    main:
    vep_annotation(vcf,params.vep_resource_dir)
    
    emit:
    vep_annotation.out.vcf
    
}


workflow qc_pretrimming {
    take:
    fasta_channel
    
    main:
    // fastqc_pre trimming
    fastqc(fasta_channel,"pre")
    //multiqc
    multiqc(fastqc.out.fastqc.collect(),"pre")
    
}

workflow trimming {
    
    take:
    fasta_channel
    
    main:
    // trim reads
    trim_illumina(fasta_channel,params.trimmomatic_file)
    
    // fastqc_post trimming
    fastqc(trim_illumina.out.trimmed_reads,"post")

    //multiqc post trimming 
    multiqc(fastqc.out.fastqc.collect(),"post")


    emit:
    trim_illumina.out.trimmed_reads
}



workflow align_and_calibrate{

    take:
    fasta_channel
    
    main:
    // split into readgroups, align using BWA MEM, mark duplicates, and set tags
    split_align_mark(fasta_channel,\
		     params.resource_dir)

    //split_align_mark.out.bam_bai.collect().view()
    // calculate alignment metrics
    alignment_metrics(split_align_mark.out.bam_bai,params.resource_dir)

    // base recalibrator
    base_recalibrator(split_align_mark.out.bam_bai,\
		      params.interval_list,\
		      params.resource_dir)

    emit:
    base_recalibrator.out.bqsr_bam_bai
}

workflow call_gvcf {

    take: bam_channel

    main:
    haplotype_caller_gvcf(bam_channel,params.interval_list,params.resource_dir)
        
    emit:
    haplotype_caller_gvcf.out.gvcf
    
}


workflow jointcall {

    take:
    gvcf

    main:
    // create channels with different chrmosomes
    chr_channel = Channel.fromList(params.joint_chromosomes )
    
    genomicsDBImport(chr_channel,gvcf,params.resource_dir)

    genomicsDBImport.out.chr_vcf.collect().view()

    // gather chromosomes
    gather(genomicsDBImport.out.chr_vcf.collect())
        
    variantRecalibratorSNP(gather.out.all_vcf_sites,\
			   params.resource_dir)
    
    variantRecalibratorINDEL(gather.out.all_vcf_sites,\
			     params.resource_dir)

    // uses as input allCHR
    applyVqsrINDEL(gather.out.all_vcf,\
		   variantRecalibratorINDEL.out.allChr_indel_recal,\
		   variantRecalibratorINDEL.out.allChr_indel_tranches)

    // uses as input indel vcf
    applyVqsrSNP(applyVqsrINDEL.out.allChr_indel_vcf,\
		 variantRecalibratorSNP.out.allChr_snp_recal,\
		 variantRecalibratorSNP.out.allChr_snp_tranches)

    emit:
    allChr_indel_vcf = applyVqsrINDEL.out.allChr_indel_vcf
    allChr_snp_vcf = applyVqsrSNP.out.allChr_snp_vcf

}
// 


workflow{

    
    // read file with fasta names
    fasta_files =  file(params.fasta).readLines()

    // create channel 
    fasta_channel = Channel.fromFilePairs(fasta_files)


    //fasta_channel.view()
    qc_pretrimming(fasta_channel)
    trimming(fasta_channel)
    
    align_and_calibrate(trimming.out)

    call_gvcf(align_and_calibrate.out)

    // collect gvcf and make joint call
    call_gvcf.out.view()
    jointcall(call_gvcf.collect())

    // annotate vcfs
    annotate(jointcall.out.allChr_indel_vcf)
    annotate(jointcall.out.allChr_snp_vcf)

    
}


