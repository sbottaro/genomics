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

include {vep_annotation } from './modules/annotation.nf'

workflow {

    //haplotype_caller_gvcf(bam_channel,params.interval_list,params.resource_dir)
    
    gvcf = Channel.fromPath(file(params.gvcfs).readLines())
    vep_annotation(gvcf,params.vep_resource_dir)

}
