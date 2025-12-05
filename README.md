# MetaGpipeline
# Metagenomics Analysis Pipeline
A comprehensive end-to-end bioinformatics pipeline for metagenomic data analysis, including quality control, assembly, binning, taxonomic profiling, functional annotation, and MAG (Metagenome-Assembled Genomes) characterization.



## Table of Contents
- [Pipeline Overview](#pipeline-overview)
- [Pipeline Modules](#pipeline-modules)
- [Output Files](#output-files)

## Pipeline Overview
This pipeline processes raw metagenomic sequencing data (paired-end FASTQ) through core steps: quality control → assembly → contig filtering → ORF prediction → MAG binning & integration → taxonomic profiling → functional annotation → gene expression quantification → MAG characterization → metabolic potential analysis.

## Pipeline Modules
### 1. Quality Control (TrimGalore)
Removes low-quality reads, adaptors, and short reads with parameters: Q30 quality threshold, minimum read length 100bp, maximum 5 N bases. Runs FastQC for quality assessment and outputs cleaned reads to `02Cleandata/Trimgalore_out/`.

### 2. Metagenomic Assembly (MegaHit)
Performs de novo assembly of metagenomic reads (minimum contig length 500bp, 30 threads). Outputs assembled contigs to `03Assembly/megahit/` and cleans up intermediate files to save space.

### 3. Contig Filtering (SeqKit)
Retains only contigs ≥1000bp for downstream analysis, outputting filtered contigs to `03Assembly/contig_1k/`.

### 4. ORF Prediction (Prodigal)
Predicts protein-coding genes from filtered contigs in metagenomic mode, generating amino acid (faa), nucleotide (fna) sequences and GFF annotations in `04ORF/prodigal/`.

### 5. MAG Binning
- **MetaWRAP**: Uses MetaBAT2, MaxBin2, and CONCOCT for multi-tool binning (30 threads, 100GB memory limit), outputting bins to `05MAG/binning/${sample}_matb2_maxb2_conc/`.
- **VAMB**: Implements variational autoencoder-based binning (sensitive to low-abundance MAGs) via contig concatenation, Minimap2 alignment, Samtools filtering, and VAMB binning, outputting to `05MAG/binning/${sample}_vamb/`.

### 6. MAG Integration & Dereplication
- **DAS Tool**: Integrates bins from multiple methods, outputting integrated bins to `05MAG/DAS/DAS_out/`.
- **dRep**: Dereplicates MAGs at 95% ANI (species-level) with filters: ≥50% completeness, ≤10% contamination. Outputs dereplicated MAGs to `05MAG/MAG_data/MAG_95/`.

### 7. Taxonomic Profiling
- **Kraken2 + Bracken**: Read-based species-level taxonomic classification with abundance correction, outputting results to `02Cleandata/reads_analysis/kraken/`.
- **mOTUs**: Marker gene-based microbial abundance profiling, outputting merged abundance tables to `02Cleandata/reads_analysis/motu/`.
- **SingleM**: 16S rRNA-based profiling at phylum/class/genus levels, outputting to `02Cleandata/reads_analysis/singlem/`.

### 8. Functional Annotation
- **KOFAMscan**: Annotates ORFs with KEGG Orthologs (E-value ≤1e-5), integrating TPM values in `04ORF/Annotation/KEGG/`.
- **MetaCerberus**: Performs multi-database annotation (KOFam_all, COG, VOG, PHROG, CAZy) of MAGs, outputting to `05MAG/MAG_data/Annotation/metacerberus/`.

### 9. Gene Expression Quantification
- **featureCounts**: Counts reads mapped to ORFs (input: BAM files + Prodigal GFF), outputting raw counts to `04ORF/Counts/`.
- **TPM Calculation**: Normalizes counts to Transcripts Per Million, outputting to `04ORF/TPM/`.

### 10. MAG Characterization
- **CoverM**: Calculates MAG coverage (TPM) from cleaned reads and dereplicated MAGs, outputting to `05MAG/MAG_data/Coverm/`.
- **GTDB-Tk**: Classifies MAGs with GTDB taxonomy and infers phylogenetic trees for Archaea/Bacteria, outputting to `05MAG/MAG_data/classify_wf_SGB/`.

### 11. MAG Metabolic Potential Analysis (METABOLIC-G)
Predicts metabolic pathways of MAGs using protein sequences from Prodigal.
- **Input**: MAG protein sequences in `05MAG/MAG_data/prodigal/`
- **Output**: Metabolic pathway analysis results in `05MAG/MAG_data/Metabolic_re/`

## Output Files
Key final outputs for downstream analysis:
| File Path | Description |
|-----------|-------------|
| `02Cleandata/reads_analysis/All_species.tsv` | Combined species-level abundance (Kraken2/Bracken) |
| `02Cleandata/reads_analysis/motu/Abundance_merge.tsv` | Merged mOTUs abundance table |
| `02Cleandata/reads_analysis/Singlem.genus_by_sample.tsv` | SingleM genus-level abundance |
| `04ORF/Annotation/KEGG/Allsample_KEGG_abundance.tsv` | Combined KO abundance + TPM |
| `05MAG/MAG_data/MAG_95/dereplicated_genomes/` | Dereplicated high/medium quality MAGs |
| `05MAG/MAG_data/classify_wf_SGB/gtdbtk.bac120.summary.tsv` | GTDB-Tk classification results |
| `05MAG/MAG_data/Annotation/metacerberus/` | Multi-database functional annotations of MAGs |
| `05MAG/MAG_data/Metabolic_re/` | MAG metabolic pathway prediction results (METABOLIC-G) |
