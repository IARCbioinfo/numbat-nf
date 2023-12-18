#!/usr/bin/env nextflow

params.label = "run_default"
params.samples = "run_default"
params.bams = 'default.bam' //To Modify, and don't forget to put index (.bai) in same directory
params.barcodes =  "../../../filtered_gene_bc_matrices/hg19/barcodes.tsv" 
params.gmap = "../../../genetic_map_hg38_withX.txt.gz" 
params.snpvcs = "../../../genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf" 
params.paneldir = "../../../1000G_hg38"
params.outdir = "res_default"
params.ncores = "1"
params.matrix = "../../../filtered_gene_bc_matrices/hg19/matrix.mtx" //To modify
params.help = null
params.eagle = "/home/claire/Eagle_v2.4.1" //To modify

process step1{

    output:
        path "${sample}_allele_counts.tsv.gz"
    publishDir "${params.outdir}", mode: "copy"

    script:
    if (params.help) {
    """
    Rscript ../../../bin/pileup_and_phase.R --help 
    """//FIND A WAY TO PRINT
    }
    else {

    """
    Rscript ../../../bin/pileup_and_phase.R --label $params.label --samples $params.samples --bams $params.bams --barcodes $params.barcodes --gmap $params.gmap --snpvcf $params.snpvcs --paneldir $params.paneldir --outdir $params.outdir --ncores $params.ncores --eagle $params.eagle
    """
    }
}

process step2{
    input:
        path step1_output

    output:
        path "*"
    publishDir "${params.outdir}", mode: "copy"
    
    script:
    //"""
    //cp ${step1_output} that_work.tsv.gz
    //"""
    """
    Rscript ../../../bin/nextflowprojet.R ${params.matrix} ${step1_output} ${params.outdir}
    """
}
workflow{
    output_step1_ch = step1()
    step2(output_step1_ch)
    }
