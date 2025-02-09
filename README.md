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

### **Step 2: Quality Control**
```bash
cd datasets/fastq_files
./run_fastqc.sh


### **Step 3: Trim Adapters**
'''bash
cd datasets/fastq_files
./run_trimming.sh

### **Step 4: Index the Genome**
'''bash
cd genome_gtf
./generate_star_index.sh

### **Step 5: Align Reads**
'''bash
cd alignment
./align_samples.sh

### **Step 6: Count Gene Expression**
'''bash
cd datasets/alignment
./feature_counts.sh


📊 Expected Output
Step	Output File	Description
FastQC	*_fastqc.html	Quality control report
Trimming	*_trimmed.fastq.gz	Cleaned reads
Alignment	*.bam	Mapped reads
featureCounts	gene_counts_matrix.txt	Gene count matrix
💡 Troubleshooting
Error Message	Possible Cause	Solution
Command not found	Missing dependency	Install the required tool (e.g., conda install ...)
Permission denied	Script lacks execution rights	Run chmod +x script.sh
Low alignment rate	Poor read quality	Check FastQC results and consider re-trimming





