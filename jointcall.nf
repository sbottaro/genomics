#!/usr/bin/env nextflow
nextflow.enable.dsl=2

VERSION="0.0"

log.info "===================================================================="
log.info " Joint call                                         "
log.info "===================================================================="

params.help = ""
if (params.help) {
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "nextflow run jointcall.nf --gvcfs list.txt "
  log.info " "
    log.info " list.txt contains gvcfs files in the following format"
    log.info " FULL_PATH/sample_01.g.vcf.gz* "
    log.info " FULL_PATH/sample_02.g.vcf.gz* "
    log.info " FULL_PATH/sample_03.g.vcf.gz* "
    log.info " "
    log.info " remember to include the * at the end so that index files will be recognized from the same path"
    log.info " "
  log.info "===================================================================="
  exit 1
}

include {haplotype_caller_gvcf; genomicsDBImport } from './modules/gatk.nf'
include {gather; variantRecalibratorSNP; variantRecalibratorINDEL } from './modules/gatk.nf'
include {applyVqsrINDEL; applyVqsrSNP } from './modules/gatk.nf'


workflow jointcall {

    take:
    gvcf

    main:
    // create channels with different chrmosomes
    chr_channel = Channel.fromPath(params.jointcall_lists )
    
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

//  The default workflow
workflow {

    //haplotype_caller_gvcf(bam_channel,params.interval_list,params.resource_dir)
    
    gvcf = Channel.fromPath(file(params.gvcfs).readLines()).collect()

    ////HaplotypeCaller
    //gvcf = Channel.fromPath(params.gvcfs).collect()
    //gvcf =  haplotype_caller_gvcf.out.gvcf.collect()
    //gvcf.collect().view()
    
    jointcall(gvcf.collect())

    //jointcall.out.view()
    jointcall.out.allChr_indel_vcf.view()
    jointcall.out.allChr_snp_vcf.view()
}
