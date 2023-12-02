### **Bulk_RNA_Seq**

# Bulk RNA-seq Data Analysis Tutorial

In this tutorial, we will guide you through the process of bulk RNA-seq data analysis (transcriptomic data analysis) using pair end samples from the NCBI GEO dataset. We will utilize the Ensembl database (GRCh38.p14) for the reference genome and reference annotation files (GFF, GTF).

For this Tutorial we have using this in the Anaconda Environment in linux

# Set up the Anaconda Environment 
1. Download [Anaconda](https://www.anaconda.com/download#downloads) from the website based on your system version
2. Installation
   ```
   bash anaconda.sh
   ```
3. Then Create a new environment:
   ```
   conda create -n Rnaseq_analysis
   ```

4. Activate the created environment
   ```
   conda activate Rnaseq_analysis
   ```
 
## Prerequisites

Make sure you have the following tools installed in the working conda environment(**Rnaseq_analysis**) :

```
conda install -c bioconda fastqc
conda install -c bioconda trimmomatic
conda install -c bioconda bwa
conda install -c bioconda subread
```

Steps to be followe




