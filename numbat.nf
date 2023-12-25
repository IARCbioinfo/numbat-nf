#!/usr/bin/env nextflow
// Copyright (C) 2023 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

nextflow.enable.dsl=2

log.info ""
log.info "--------------------------------------------------------"
log.info "  numbat-nf v1.0: CNV calling from single-cell or spatial RNA-seq"
log.info "--------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

params.label = "run_default"
params.samples = "run_default"
params.input_file = null
params.bams = null //'default.bam' //To Modify, and don't forget to put index (.bai) in same directory
params.barcodes = null //"cell_barcodes_v2.txt" 
params.gmap = "/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz" 
params.snpvcf = "/data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf" 
params.paneldir = "/data/1000G_hg38"
params.output_folder = "numbat_results"
params.cpu = "8"
params.mem = 10
params.matrix = null //"filtered_gene_bc_matrices/hg38/matrix.mtx" //To modify
params.help = null
params.eagle = "/home/Eagle_v2.4.1" //To modify


if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run iarcbioinfo/numbat-nf [-with-docker] --input_file input.tsv [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--bams                      <PATH>                      <DESCRIPTION>"
    log.info ""
    log.info "Optional arguments:"
    log.info '--output_folder  STRING                 Output folder (default: .).'
    log.info '--cpu            INTEGER                Number of cpu used by bwa mem and sambamba (default: 8).'
    log.info '--mem            INTEGER                Size of memory used for mapping (in GB) (default: 10).'
    log.info ""
    log.info "Flags:"
    log.info "--<FLAG>                                                    <DESCRIPTION>"
    log.info ""
    exit 0
} else {
/* Software information */
log.info "help:                               ${params.help}"
log.info "bams   = ${params.bams}"
log.info "output_folder    = ${params.output_folder}"
log.info "barcodes         = ${params.barcodes}"
log.info "matrix           = ${params.matrix}"
}

process pileup_and_phase{
    cpus params.cpu
    memory params.mem+'GB'    
    tag {ID}

    input:
    tuple val(ID), path(matrix), path(barcodes), path(bam), path(bai)
    path gmap
    path snpvcf
    path paneldir

    output:
        tuple path("*_allele_counts.tsv.gz"), path(matrix)
        tuple path("phasing*"), path("pileup")
    publishDir "${params.output_folder}/intermediate/allele_counts", mode: "copy", pattern: "*_allele_counts.tsv.gz"
    publishDir "${params.output_folder}/intermediate", mode: "copy", pattern: "phasing*"

    script:
    """
    Rscript $projectDir/bin/pileup_and_phase.R --label $ID --samples $ID --bams $bam --barcodes $barcodes --gmap $gmap --snpvcf $snpvcf --paneldir $paneldir --outdir "./" --ncores $params.cpu --eagle $params.eagle
    mv phasing.log phasing/${ID}_phasing.log
    """
}

process numbat{
    cpus params.cpu
    memory params.mem+'GB'

    input:
        tuple path(allele_counts), path(matrix)

    output:
        path "*"
    publishDir "${params.output_folder}/results", mode: "copy"
    
    script:
    """
    Rscript $projectDir/bin/nextflowprojet.R $matrix $allele_counts "./" $params.cpu
    """
}

workflow{
    if(params.input_file){
	bams = Channel.fromPath("${params.input_file}")
     	          .splitCsv( header: true, sep: '\t', strip: true )
	       .map { row -> [ row.ID , row.matrix , row.barcodes, file(row.bam), file(row.bam+'.bai') ] }
	       .view()
    }else{
        bams     = Channel.fromFilePairs( params.bams+'{,.bai}' ).view() //to modify if needed to concatenate barcondes and matrix
        matrix   = Channel.fromPath(params.matrix)
        barcodes = Channel.fromPath(params.barcodes)
    }
    gmap     = file(params.gmap)
    snpvcf   = file(params.snpvcf)
    paneldir = file(params.paneldir)

    output_step1_ch = pileup_and_phase(bams,gmap,snpvcf,paneldir)
    numbat(output_step1_ch[0])
    }
