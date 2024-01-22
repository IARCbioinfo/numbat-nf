# numbat-nf
## A bioinformatic workflow to call copy number variants reconstruct cell phylogenies from single-cell or spatial transcriptomics data
## Description
This is a workflow in nextflow to preprocess and use the numbat software for the analysis of single-cell or spatial transcromics data to study Copy number variations (CNVs).

## Dependencies

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.
2. External software:
- [numbat](https://kharchenkolab.github.io/numbat/articles/numbat.html)

You can avoid installing all the external software by only installing Docker, and downloading the [numbat docker file](https://kharchenkolab.github.io/numbat/articles/numbat.html). See the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository for more information.

## Input
  | Type      | Description     |
  |-----------|---------------|
  | input_file    | Tab-separated file with columns ID (sample ID), matrix_folder (cellranger output folder with expression matrix .mtx, barcode, and features files), and bam (alignment files) |

## Parameters

  * #### Mandatory
| Name      | Example value | Description     |
|-----------|---------------|-----------------|
|--gmap | /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz |          Path to genetic map provided by Eagle2 |
|--eagle | eagle |        Path to Eagle2 binary file  |
|--snpvcf | /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf |      SNP VCF for pileup (e.g. genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf) |
|--paneldir | /data/1000G_hg38 |  Directory to phasing reference panel (e.g. 1000G_hg38) |

All these files can be downloaded from the numbat website, or are contained in the numbat docker container (note: the default paths correspond to the location within the docker container).

  * #### Optional
| Name      | Default value | Description     |
|-----------|---------------|-----------------|
|--ncores | 8 |      Number of cores to use.|
|--mem | 10 |      Memmory (in Gb) to use.|

## Usage
```nextflow run script_nf.nf --input_file input.tsv [--OPTIONS] OPTION```

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------|
  | Nicolas Alcala*    | alcalan@iarc.who.int | Developer to contact for support |
  | Yanis Sindt-Baret | Developer |
  | Quentin Ohayon    |  | Developer |
  | Natacha Doutrelea |  | Developer |
  | Claire Berthaud   |  | Developer |
  
## FAQ
### Prerequisites
To work, the bam file needs to come with its index file (.bai); if necessary, it can be generated using the command:
```samtools index -b BAME_FILE.bam```
