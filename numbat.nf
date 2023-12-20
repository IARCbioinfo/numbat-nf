#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.label = "run_default"
params.samples = "run_default"
params.bams = 'default.bam' //To Modify, and don't forget to put index (.bai) in same directory
params.barcodes =  "cell_barcodes_v2.txt" 
params.gmap = "/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz" 
params.snpvcf = "/data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf" 
params.paneldir = "/data/1000G_hg38"
params.output_folder = "numbat_results"
params.cpu = "8"
params.matrix = "filtered_gene_bc_matrices/hg38/matrix.mtx" //To modify
params.help = null
params.eagle = "/home/Eagle_v2.4.1" //To modify


log.info "bams   = ${params.bams}"
log.info "output_folder    = ${params.output_folder}"
log.info "barcodes         = ${params.barcodes}"
log.info "matrix           = ${params.matrix}"

process step1{
    cpus params.cpu

    input:
    tuple val(ID), path(bambai)
    path matrix
    path barcodes
    path gmap
    path snpvcf
    path paneldir

    output:
        path "${sample}_allele_counts.tsv.gz"
    publishDir "${params.output_folder}", mode: "copy"

    script:
    if (params.help) {
    """
    Rscript $projectDir/bin/pileup_and_phase.R --help 
    """//FIND A WAY TO PRINT
    }
    else {
    """
    Rscript $projectDir/bin/pileup_and_phase.R --label $params.label --samples $params.samples --bams ${bambai[0]} --barcodes $barcodes --gmap $gmap --snpvcf $snpvcf --paneldir $paneldir --outdir "./" --ncores $params.cpu --eagle $params.eagle
    """
    }
}

process step2{
    cpus '2'

    input:
        path step1_output

    output:
        path "*"
    publishDir "${params.output_folder}", mode: "copy"
    
    script:
    """
    Rscript $projectDir/bin/nextflowprojet.R ${params.matrix} ${step1_output} "./"
    """
}

workflow{
    bams     = Channel.fromFilePairs( params.bams+'{,.bai}' ).view()
    matrix   = Channel.fromPath(params.matrix)
    barcodes = Channel.fromPath(params.barcodes)
    gmap     = file(params.gmap)
    snpvcf   = file(params.snpvcf)
    paneldir = file(params.paneldir)

    output_step1_ch = step1(bams,matrix,barcodes,gmap,snpvcf,paneldir)
    step2(output_step1_ch)
    }
