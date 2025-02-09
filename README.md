# RNA-Seq Pipeline for LBD Analysis

This repository contains a **bioinformatics pipeline for processing RNA-Seq data**, from raw reads to **differential gene expression (DEG) analysis**. It covers quality control, read trimming, alignment, gene counting, and enrichment analysis.

---

## 📌 Table of Contents
- [📜 Overview](#-overview)
- [🛠 Requirements](#-requirements)
- [📂 Directory Structure](#-directory-structure)
- [⚡ Running the Pipeline](#-running-the-pipeline)
  - [Step 1: Download Raw Data](#step-1-download-raw-data)
  - [Step 2: Quality Control](#step-2-quality-control)
  - [Step 3: Trim Adapters](#step-3-trim-adapters)
  - [Step 4: Index the Genome](#step-4-index-the-genome)
  - [Step 5: Align Reads](#step-5-align-reads)
  - [Step 6: Count Gene Expression](#step-6-count-gene-expression)
- [📊 Expected Output](#-expected-output)
- [💡 Troubleshooting](#-troubleshooting)
- [📜 Citation](#-citation)

---

## 📜 Overview
This pipeline is designed to analyze **RNA-Seq data related to LBD (Lewy Body Dementia)**. The workflow consists of:
- **Quality control** using FastQC
- **Adapter trimming** using Trimmomatic
- **Alignment** using STAR
- **Gene expression quantification** using featureCounts
- **Downstream analysis** with R (DESeq2, GO enrichment)

---

## 🛠 Requirements
### **Software**
Make sure the following dependencies are installed:

| Tool         | Purpose                  | Installation Command                      |
|--------------|--------------------------|-------------------------------------------|
| FastQC       | Quality control          | `sudo apt install fastqc`                |
| Trimmomatic  | Adapter trimming         | `wget https://github.com/usadellab/Trimmomatic` |
| STAR         | Read alignment           | `conda install -c bioconda star`         |
| samtools     | BAM processing           | `conda install -c bioconda samtools`     |
| featureCounts| Gene counting            | `conda install -c bioconda subread`      |
| R + DESeq2   | Differential Expression  | `R -e 'install.packages("DESeq2")'`      |

---

## 📂 Directory Structure
![Directory Structure](https://github.com/Djinho/rna-seqLDB/blob/main/directory_structure.png)

---

## ⚡ Running the Pipeline

### **Step 1: Download Raw Data**
```bash
cd datasets
bash SRA_2_FAST.sh





Step 4: Index the Genome
cd genome_gtf
chmod +x generate_star_index.sh
./generate_star_index.sh
****











