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

| Tool         | Purpose                  | Installation Command |
|-------------|--------------------------|----------------------|
| FastQC      | Quality control           | `sudo apt install fastqc` |
| Trimmomatic | Adapter trimming          | `wget https://github.com/usadellab/Trimmomatic` |
| STAR        | Read alignment            | `conda install -c bioconda star` |
| samtools    | BAM processing            | `conda install -c bioconda samtools` |
| featureCounts | Gene counting           | `conda install -c bioconda subread` |
| R + DESeq2  | Differential Expression   | `R -e 'install.packages("DESeq2")'` |

---

## 📂 Directory Structure

rna-seqLDB/
├── 📂 alignment/ # STAR alignments
│   ├── 📄 align_samples.sh
│   ├── 📄 feature_counts.sh
│   ├── 📄 gene_counts_matrix.txt
│   ├── 📄 *.bam (Aligned BAM Files)
│   └── 📄 *.bai (BAM Index Files)
├── 📂 datasets/ # Data & Scripts
│   ├── 📄 SRA_2_FAST.sh
│   ├── 📂 fastq_files/ # Raw FASTQ Data
│   │   └── 📄 *.fastq.gz
│   ├── 📂 fastqc_results/ # FastQC Outputs
│   │   ├── 📄 *.html (Reports)
│   │   └── 📄 *.zip (Raw Data)
│   ├── 📂 trimmed/ # Trimmed Reads
│   │   ├── 📄 *_paired.fastq.gz
│   │   └── 📄 *_unpaired.fastq.gz
│   ├── 📄 run_fastqc.sh
│   └── 📄 run_trimming.sh
├── 📂 genome_gtf/ # Genome Files
│   ├── 📄 Homo_sapiens.GRCh38.113.gtf
│   ├── 📄 Homo_sapiens.GRCh38.dna.primary_assembly.fa
│   └── 📂 STAR_index/
│       ├── 📄 Genome
│       ├── 📄 SA
│       ├── 📄 chrName.txt
│       └── 📄 sjdbList.out.tab
└── 📂 scripts/ # Custom Scripts
    ├── 📄 generate_star_index.sh
    └── 📄 investigate_env.sh

---

## ⚡ Running the Pipeline

### **Step 1: Download Raw Data**
```bash
cd datasets
bash SRA_2_FAST.sh
📝 What it does: Downloads sequencing data from SraAccList.csv.

Step 2: Quality Control
bash
Copy
Edit
cd datasets/fastq_files
bash run_fastqc.sh
📝 What it does: Runs FastQC and saves results in fastqc_results/.

Step 3: Trim Adapters
bash
Copy
Edit
cd datasets/fastq_files
bash run_trimming.sh
📝 What it does: Removes adapter sequences using Trimmomatic.

Step 4: Index the Genome
bash
Copy
Edit
cd genome_gtf
bash generate_star_index.sh
📝 What it does: Generates a STAR genome index (needed for alignment).

Step 5: Align Reads
bash
Copy
Edit
cd alignment
bash align_samples.sh
📝 What it does: Aligns reads to the GRCh38 reference genome.

Step 6: Count Gene Expression
bash
Copy
Edit
cd datasets/alignment
bash feature_counts.sh
📝 What it does: Runs featureCounts to generate a gene count matrix.

📊 Expected Output
Step	Output File	Description
FastQC	*_fastqc.html	Quality control report
Trimming	*_trimmed.fastq.gz	Cleaned reads
Alignment	*.bam	Mapped reads
featureCounts	gene_counts_matrix.txt	Gene count matrix
💡 Troubleshooting
🔹 Common Issues & Fixes
Error Message	Possible Cause	Solution
Command not found	Missing dependency	Install the required tool (conda install ...)
Permission denied	Script lacks execution rights	Run chmod +x script.sh
Low alignment rate	Poor read quality	Check FastQC results and consider re-trimming
📜 Citation
If you use this pipeline, please cite:

STAR: Dobin et al., 2013
featureCounts: Liao et al., 2014
Trimmomatic: Bolger et al., 2014
📌 Additional Notes
Modify file paths before running the scripts.
Adjust thread count (-T 4) based on your system.





