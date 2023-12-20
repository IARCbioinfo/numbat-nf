# PROJECT 4: Implement a bioinformatic workflow to call copy number variants reconstruct cell phylogenies from single-cell sequencing data

Goal of this project is to implement a workflow in nextflow to preprocess and use numbat software for the analysis of single-cell sequencing data to study Copy number variations (CNVs).

To use this workflow, please follow these steps. 

## Step 0: Installation

Follow the step **Installation** on this page :
https://kharchenkolab.github.io/numbat/articles/numbat.html

## Step 1: Build bam index (if needed)
To work, the bam need to come with his index file (.bai), to do that you can use this command:

```samtools index -b BAME_FILE.bam```

Please change BAME_FILE by the file you wanna use and make sure you have samtools.

## Step 2: Launch workflow
```nextflow run script_nf.nf [--OPTIONS] OPTION```

Using the following documentary, all parameter are necessary.
```
--label LABEL        Individual label. One per run.
--samples SAMPLES    Sample name
--bams BAMS          BAM file
--barcodes BARCODES  Cell barcode file
                       (e.g. "filtered_gene_bc_matrices/hg19/barcodes.tsv")
--gmap GMAP          Path to genetic map provided by Eagle2 (e.g.
                       Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz)
--eagle EAGLE        Path to Eagle2 binary file 
--snpvcf SNPVCF      SNP VCF for pileup (e.g. genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf)
--paneldir PANELDIR  Directory to phasing reference panel (e.g. 1000G_hg38)
--outdir OUTDIR      Output directory
--matrix MATRIX      Gene x cell integer UMI count matrix
                       (e.g. filtered_gene_bc_matrices/hg19/matrix.mtx)
--ncores NCORES      Number of cores to use.
```

You may change in script_nf.nf the default value for :
  - params.gmap
  - params.snpvcs
  - params.paneldir
  - params.ncores
  - params.eagle
    
To your locals path, then it can be fixed for any run with these references. 

For example, if your local path are fixe, the bash command may look like : 

```nextflow run script_nf.nf --samples "example" --bams "example.bam" --barcodes "barcodes.bam" --matrix "matrix.bam" --outdir "result"```

## Step 3: Results
