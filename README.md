## Introduction
This repository contains data and code for our analysis of the genome of a Plasmodium falciparum strain with reduced susceptibility to Artemisinin-Based Combination Therapy (ACT) drugs. We performed whole-genome sequencing of the HSOG3 clinical isolate and compared it with a large dataset of global P. falciparum samples. We identified and annotated SNPs that are associated with drug resistance and phylogeography, and calculated the allele frequencies for each population.

## Methods
### Analysis of NGS data
We used a variant calling pipeline for P. falciparum genome sequencing data, adapted from the gencorefacility/variant-calling-pipeline-gatk4 GitHub repository. The pipeline was implemented in Nextflow version 21.10.6 and run with Docker to ensure reproducibility. The pipeline consists of several steps: 
- Quality control of raw reads using FastQC and MultiQC
- Alignment of reads to the reference genome (Pf3D7_all_v3.fasta) using BWA-MEM
- Marking of PCR duplicates using Picard MarkDuplicates
- Base quality score recalibration using GATK BaseRecalibrator
- Variant calling using GATK HaplotypeCaller version 4.1.9.0
- Variant filtering using GATK VariantFiltration and custom filters

We created the index files for the reference genome (BWA index, SAMtools index, and Picard dictionary) using the broadinstitute/gatk Docker image. We modified the nextflow.config file to specify the input files, output directory, and other parameters. We ran the pipeline using the -with-docker option to use the gencorefacility/variant-calling-pipeline-gatk4 Docker image. The resulting VCF file was filtered using several quality control filters including FS_filter, MQRankSum_filter, MQ_filter, QD_filter, ReadPosRankSum_filter, and SOR_filter to remove low-quality variants. SnpEff version 4.3i was used to annotate the VCF file.

### Annotating SNPs related to drug resistance and phylogeography
We used custom Python code (version 3.10.12) to analyze the SNPs present in the P. falciparum genome sequencing data. For drug resistance-related SNPs, we manually curated a list of genes from the literature review. For phylogeographic informative SNPs, we used the list from malaria-db GitHub repository, which contains SNPs that indicate the geographic origin of the parasite strains. We used scikit-allel (version 1.3.7) for efficient reading of the VCF file containing our sample's variant calls. We also used pandas (version 2.1.1) and numpy (version 1.25.2) for data manipulation and visualization.

### Analysis of the "Pf7 Dataset"
We used the "Pf7 dataset" from the MalariaGEN P. falciparum Community Project, which can be accessed at ftp://ngs.sanger.ac.uk/production/malaria/Resource/34/Pf7.zarr.zip or from the malariagen_data python library. This dataset contains genotype calls on 10,145,661 SNPs and short indels across 20,864 samples from 33 countries, structured in an xarray format. Out of the 16,203 whole-genome sequencing samples that passed MalariaGEN's quality control measures, we selected a subset of 10,348 samples that had an Fws value of ≥ 0.95. The Fws values for all samples were retrieved from https://www.malariagen.net/sites/default/files/Pf7_fws.txt. We used custom Python code to select the samples and specific coordinates of interest (chromosome and genomic position). The population annotation for each sample was obtained from https://www.malariagen.net/sites/default/files/Pf7_samples.txt. We calculated the percentage of the alternative allele for each population.

## Data
- **GATK_ins**: This directory houses the input files required for the Genome Analysis Toolkit (GATK), a software package designed for the analysis of high-throughput sequencing data. It includes various formats of the reference genome (Pf3D7_all_v3), along with scripts and configuration files.

- **GATK_outs**: This directory stores the output files generated by the GATK analysis pipeline. It contains the file 'hsog3catall_filtered_snps.ann.vcf', which is a Variant Call Format (VCF) file comprising SNP calls from the HSOG3 isolate after filtering and annotation. Please note that to ensure patient confidentiality and privacy, the raw sequencing data obtained from the clinical isolate cannot be shared publicly. It also contains the file 'magical_franklin_report.csv', which is a report generated by the GATK pipeline that provides a comprehensive overview of the sequencing data and variant calling results.

- **HSOG_Pf7_SNP_analysis**: This directory encompasses files related to the SNP analysis of the HSOG clinical isolate. The 'DR_genes.csv' file provides information about the list of drug resistance genes curated for this study. The 'geo_barcode.tsv' file, obtained from the malaria-profiler database, contains information about geographically informative SNPs. The 'Pf7_sample_index_FWS_095.txt' file includes sample indexing information for the Pf7 dataset samples with an Fws value of ≥ 0.95. The Python scripts 'pf7.py' and 'vcf_subset.py' contain code used for data analysis.

## How to reference
If you use this in your project, please cite it as follows:

BibTex format

```bibtex
@article{Casanova_2023, 
    title={Artemisinin resistance-associated gene mutations in Plasmodium falciparum: A case study of severe malaria from Mozambique}, 
    ISSN={1477-8939}, 
    url={http://dx.doi.org/10.1016/j.tmaid.2023.102684}, 
    DOI={10.1016/j.tmaid.2023.102684}, 
    journal={Travel Medicine and Infectious Disease}, 
    publisher={Elsevier BV}, 
    author={Casanova, Daniela and Baptista, Vitória and Costa, Magda and Freitas, Bruno and Pereira, Maria and Calçada, Carla and Mota, Paula and Kythrich, Olena and Pereira, Maria Helena Jacinto Sarmento and Osório, Nuno S. and Veiga, M Isabel}, 
    year={2023}, 
    month={dec},
}
```

APA format
```
Casanova, D., Baptista, V., Costa, M., Freitas, B., Pereira, M., ... Osório, N. S. & Veiga, M. I. (2023). Artemisinin resistance-associated gene mutations in Plasmodium falciparum: A case study of severe malaria from Mozambique. Travel Medicine and Infectious Disease, 102684. doi:10.1016/j.tmaid.2023.102684
```

[![Link for the publication](https://images.app.goo.gl/XNCyWYaEG8CVPpa76)](http://dx.doi.org/10.1016/j.tmaid.2023.102684)

## Acknowledgements
We extend our gratitude to all individuals who contributed directly or indirectly to this project, providing valuable insights and expertise.
