#!/bin/bash
# Comprehensive Metagenomics Analysis Pipeline
# Author: Zhihao Zhang
# Date: 2025.11.30
# Description: End-to-end pipeline for metagenomic data processing including:
# - Quality control and trimming
# - Assembly
# - ORF prediction
# - MAG binning (MetaBAT2, MaxBin2, CONCOCT, VAMB)
# - DAS Tool integration
# - Taxonomic profiling (Kraken2/Bracken, mOTUs, SingleM)
# - Functional annotation (KEGG, COG, CAZy, PHROG)
# - MAG quality assessment and taxonomic classification (GTDB-Tk)
# - Expression quantification (featureCounts/TPM)
# - Coverage calculation (CoverM)

###########################################################
# Step 1: Quality Control with TrimGalore
###########################################################
# Load TrimGalore environment
source /datanode02/zhangzh/.apps/trimgalore.sh 

# Process each sample with quality filtering
# --paired: Paired-end reads processing
# --quality 30: Minimum PHRED quality score
# --length 100: Minimum read length after trimming
# --max_n 5: Maximum N bases allowed
# -j 8: Threads for FastQC
# --fastqc: Run FastQC for quality assessment
# -o: Output directory
# &>: Redirect both stdout and stderr to log file
for i in  `cat sample.name`; do
    trim_galore --paired --quality 30 --length 100 --max_n 5 -j 8 --fastqc  \
        -o 02Cleandata/Trimgalore_out  \
        01Rawdata/${i}_1.fq.gz 01Rawdata/${i}_2.fq.gz &> 02Cleandata/Trimgalore_out/${i}.log
done 

# Rename trimmed files by removing "_val" suffix from filenames
# Process forward reads
\ls 02Cleandata/Trimgalore_out/*gz  |  
    awk -F"/" '$4=$3{sub("_1_val","",$3); print "mv "$1,$2,$4"\t"$1,$2,$3 }' OFS="/" | sh

# Process reverse reads
\ls 02Cleandata/Trimgalore_out/*gz  |  
    awk -F"/" '$4=$3{sub("_2_val","",$3); print "mv "$1,$2,$4"\t"$1,$2,$3 }' OFS="/" | sh

###########################################################
# Step 2: Metagenomic Assembly with MegaHit
###########################################################
# Load MegaHit environment
source /apps/source/megahit-1.2.9.sh

# Assemble each sample with MegaHit
# -t 30: Threads
# --min-contig-len 500: Minimum contig length
# -o: Output directory
# --out-prefix: Output prefix
# -1/-2: Input paired-end reads
for i in `cat sample.name`; do
    megahit  -t 30  --min-contig-len 500 \
        -o 03Assembly/megahit/${i}   --out-prefix  ${i} \
        -1 02Cleandata/Trimgalore_out/${i}_1.fq.gz    \
        -2 02Cleandata/Trimgalore_out/${i}_2.fq.gz
done 

# Clean up intermediate files to save space
for i in `cat sample.name`; do
    rm  03Assembly/megahit/${i}/intermediate_contigs -r
done 

# Add sample name prefix to contig IDs for traceability
# Prevents duplicate contig IDs across samples
for i in 03Assembly/megahit/*/*fa; do 
    n=$(basename $i .contigs.fa)
    # Check if sample prefix already exists to avoid duplication
    if ! grep -q "^>${n}_" "$i"; then
        sed "s:>:>${n}_:g" "$i" -i
    fi
done

###########################################################
# Step 3: Filter Contigs by Length (≥1000bp)
###########################################################
# Load SeqKit environment
source /apps/source/seqkit-2.2.0.sh 

# Filter contigs to retain only those ≥1000bp
# -m 1000: Minimum sequence length
for i  in 03Assembly/megahit/*/*fa ; do
    n=$(basename $i .contigs.fa)
    seqkit seq -m 1000  $i > 03Assembly/contig_1k/${n}_contigs1k.fa
done
 
###########################################################
# Step 4: ORF Prediction with Prodigal
###########################################################
# Load Prodigal environment
source /apps/source/prodigal-2.6.3.sh 

# Predict ORFs from filtered contigs (meta mode for metagenomes)
# -p meta: Metagenomic mode (mixed genomes)
# -q: Quiet mode
# -m: Treat translation table as minimal (stop codon only)
# -i: Input contig file
# -a: Output amino acid sequences (faa)
# -d: Output nucleotide sequences (fna)
# -o: Output GFF annotation file
# -f gff: Output format (GFF)
for i in `cat sample.name` ; do
    prodigal -p meta -q -m \
        -i  03Assembly/contig_1k/${i}_contigs1k.fa   \
        -a 04ORF/prodigal/${i}_contigs1k_prodigal.faa \
        -d 04ORF/prodigal/${i}_contigs1k_prodigal.fna \
        -o 04ORF/prodigal/${i}_contigs1k_prodigal.gff \
        -f gff
done 

###########################################################
# Step 5: MAG Binning with MetaWRAP (MetaBAT2, MaxBin2, CONCOCT)
###########################################################
# Load MetaWRAP environment
source /datanode02/zhangzh/.apps/metawrap.sh

# Run multiple binning algorithms with MetaWRAP
# --metabat2 --maxbin2 --concoct: Enable all three binning methods
# -o: Output directory
# -t 30: Threads
# -m 100: Memory limit (GB)
# -a: Assembled contigs (≥1000bp)
# Last arguments: Paired-end reads for coverage calculation
for i in `cat sample.name` ; do
    metawrap binning --metabat2 --maxbin2 --concoct \
        -o 05MAG/binning/${i}_matb2_maxb2_conc  \
        -t 30 -m 100  \
        -a   03Assembly/contig_1k/${i}_contigs1k.fa    \
        02Cleandata/Trimgalore_out/${i}_1.fq.gz     \
        02Cleandata/Trimgalore_out/${i}_2.fq.gz
done 

###########################################################
# Step 6: MAG Binning with VAMB (Variational Autoencoder)
###########################################################
# Load VAMB and Samtools environments
source /datanode02/zhangzh/.apps/vamb.sh
source /apps/source/samtools-1.15.1.sh

# VAMB binning pipeline (more sensitive for low-abundance MAGs)
# Steps:
# 1. Concatenate contigs
# 2. Create minimap2 index
# 3. Align reads to contigs
# 4. Filter alignments (remove unmapped/secondary)
# 5. Run VAMB binning
for i in `cat sample.name` ; do
    concatenate.py 05MAG/binning/${i}.fa.gz  03Assembly/contig_1k/${i}_contigs1k.fa   &>05MAG/binning/${i}_concatenate.log &&  
    minimap2 -d 05MAG/binning/${i}.fa.mmi 05MAG/binning/${i}.fa.gz &> 05MAG/binning/${i}_index.log &&  
    minimap2 -t 30 -N 50 -ax sr  05MAG/binning/${i}.fa.mmi   \
        02Cleandata/Trimgalore_out/${i}_1.fq.gz     \
        02Cleandata/Trimgalore_out/${i}_2.fq.gz |    
        samtools view -F 3584 -b --threads 30 > 05MAG/binning/${i}.bam  &&   
    time vamb --outdir 05MAG/binning/${i}_vamb \
        --fasta 05MAG/binning/${i}.fa.gz \
        --bamfiles 05MAG/binning/${i}.bam \
        -o C --minfasta 200000
done 

# Clean up BAM files to save space
for i  in `cat sample.name` ; do
    rm 05MAG/binning/${i}.bam
done 

###########################################################
# Step 7: Rename MAG Files for Consistency
###########################################################
cd 05MAG

# Add sample-specific prefixes to MAG filenames from different binning methods
# _ma_: MaxBin2
# _me_: MetaBAT2
# _co_: CONCOCT
for i in  binning/*/maxbin2_bins/* ; do
    n=${i#*/}
    mv  ${i}  ${i%/*}/${n//_matb2_maxb2_conc\/maxbin2_bins\//_ma_}
done 

for i in  binning/*/metabat2_bins/* ; do
    n=${i#*/}
    mv  ${i}  ${i%/*}/${n//_matb2_maxb2_conc\/metabat2_bins\//_me_}
done

for i in  binning/*/concoct_bins/* ; do
    n=${i#*/}
    mv  ${i}  ${i%/*}/${n//_matb2_maxb2_conc\/concoct_bins\//_co_}
done

# Rename VAMB bins
for i in  binning/*vamb/bins/* ; do
    n=${i#*/}
    mv ${i}  ${i%/*}/${n//_vamb\/bins\//_}
done 

# Clean up VAMB contig ID prefixes
for i in  binning/*vamb/bins/* ; do
    sed 's/>S1C/>/g' ${i} -i
done 

cd .. 

###########################################################
# Step 8: Integrate Bins with DAS Tool
###########################################################
# Load DAS Tool environment
source /datanode02/zhangzh/.apps/DAS_tool.sh

# Create contig-to-bin mapping files for each binning method
# Required input format for DAS Tool
for i in `cat sample.name`; do
    /datanode02/zhangzh/minisoft/DAS_Tool/src/Fasta_to_Contig2Bin.sh \
        -e fa \
        -i 05MAG/binning/${i}_matb2_maxb2_conc/maxbin2_bins/ > 05MAG/DAS/DAS_data/${i}_maxbin.tsv
done 

for i in `cat sample.name`; do
    /datanode02/zhangzh/minisoft/DAS_Tool/src/Fasta_to_Contig2Bin.sh \
        -e fa \
        -i 05MAG/binning/${i}_matb2_maxb2_conc/metabat2_bins/ > 05MAG/DAS/DAS_data/${i}_metabat.tsv
done 

for i in `cat sample.name`; do
    /datanode02/zhangzh/minisoft/DAS_Tool/src/Fasta_to_Contig2Bin.sh \
        -e fa \
        -i 05MAG/binning/${i}_matb2_maxb2_conc/concoct_bins/ > 05MAG/DAS/DAS_data/${i}_concot.tsv
done 

for i in `cat sample.name`; do
    /datanode02/zhangzh/minisoft/DAS_Tool/src/Fasta_to_Contig2Bin.sh \
        -e fna \
        -i 05MAG/binning/${i}_vamb/bins/  > 05MAG/DAS/DAS_data/${i}_vamb.tsv
done 

# Load OrthoFinder (required for DAS Tool protein clustering)
source /datanode02/zhangzh/.apps/orthofinder.sh

# Run DAS Tool to integrate bins from multiple methods
# -i: Comma-separated contig-to-bin mappings
# -l: Labels for each binning method
# --proteins: Protein sequences for quality assessment
# -c: Contig fasta file
# -o: Output directory
# --search_engine diamond: Use DIAMOND for homology search
# --write_bins: Output final integrated bins
# --score_threshold 0: No score threshold (retain all bins)
for i in `cat sample.name`; do
    DAS_Tool -i 05MAG/DAS/DAS_data/${i}_maxbin.tsv,05MAG/DAS/DAS_data/${i}_metabat.tsv,05MAG/DAS/DAS_data/${i}_concot.tsv,05MAG/DAS/DAS_data/${i}_vamb.tsv  \
        -l maxbins,metabat2,concoct,vamb \
        --proteins 04ORF/prodigal/${i}_contigs1k_prodigal.faa  \
        -c  03Assembly/contig_1k/${i}_contigs1k.fa  \
        -o 05MAG/DAS/DAS_out/${i}_output/ \
        --search_engine diamond \
        --write_bins   --score_threshold 0
done 

# Collect all DAS Tool integrated bins into one directory
mv 05MAG/DAS/DAS_out/*/_DASTool_bins/*  05MAG/MAG_data/DAS_bins/

###########################################################
# Step 9: Dereplicate MAGs with dRep
###########################################################
cd 05MAG/MAG_data

# Load CheckM environment (required for dRep quality assessment)
source  /datanode02/zhangzh/.apps/checkm.sh
export CHECKM_DATA_PATH="/datanode02/zhangzh/database/Checkm_database"

# Dereplicate MAGs at 95% ANI (species-level)
# -pa 0.95: Primary ANI threshold (95%)
# -cm larger: Keep larger MAG when duplicates found
# -p 30: Threads
# -sa 0.95: Secondary ANI threshold
# -comp 50: Minimum completeness (50%)
# -con 10: Maximum contamination (10%)
# -g: Input MAG files
time dRep dereplicate MAG_95  -pa 0.95 -cm larger -p 30  \
    -sa 0.95 -comp 50 -con 10 -g DAS_bins/*.fa

cd  ../../ 

###########################################################
# Step 10: Taxonomic Profiling (Kraken2/Bracken)
###########################################################
# Load Kraken2 environment
source /datanode02/zhangzh/.apps/kraken2.sh

# Run Kraken2 for read-based taxonomic classification
# --db: Kraken2 database path
# --threads 10: Threads
# --report: Taxonomic report file
# --report-minimizer-data: Include minimizer data in report
# --output: Raw classification output
# --paired: Paired-end reads
for i in `cat sample.name`; do
    kraken2 --db /datanode02/zhangzh/database/kraken2_database  \
        --threads 10 \
        --report 02Cleandata/reads_analysis/kraken/${i}.report  \
        --report-minimizer-data \
        --output 02Cleandata/reads_analysis/kraken/${i}.output  \
        --paired 02Cleandata/Trimgalore_out/${i}_1.fq.gz  \
        02Cleandata/Trimgalore_out/${i}_2.fq.gz
done

# Load Bracken environment (conda)
source activate /datanode02/zhangzh/.conda/envs/palmid
cd 02Cleandata/reads_analysis

# Run Bracken to adjust Kraken2 abundances (species-level)
# -d: Database path
# -i: Input Kraken2 report
# -o: Output Bracken file
# -w: Output Bracken report
# -r 150: Read length
# -l S: Taxonomic level (Species)
for i in ./kraken/*report; do
    n=$(basename $i .report) ;
    /datanode02/zhangzh/minisoft/Bracken/bracken -d /datanode02/zhangzh/database/kraken2_database  \
        -i ${i} \
        -o bracken/${n}.S.bracken \
        -w bracken/${n}.bracken.report \
        -r 150 -l S
done 

# Clean up large Kraken2 raw output files
for i in `cat sample.name`; do
    rm 02Cleandata/reads_analysis/kraken/${i}.output
done 

# Convert Bracken reports to MPA format (compatible with downstream analysis)
for i in ./kraken/*report; do
    n=$(basename $i .report) ;
    python  /datanode02/zhangzh/minisoft/KrakenTools/kreport2mpa.py \
        -r bracken/${n}.bracken.report   \
        -o bracken/${n}.tax
done

# Extract Kingdom-level taxonomic data
for i in bracken/*tax; do
    n=$(basename $i .tax);
    awk -F"|" '$NF~/k__/{print $0}' ${i} | sed "1i Kindom\t${n}" > kra_tax/${n}_kindom.tsv
done 

# Extract Species-level taxonomic data
for i in bracken/*tax; do
    n=$(basename $i .tax);
    awk -F"|" '$NF~/s__/{print $0}' ${i} | sed "1i Species\t${n}" > kra_tax/${n}_species.tsv
done

# Add sample name prefix to species data
for i in kra_tax/*_species.tsv ; do
    n=$(basename $i _species.tsv) ;
    sed "s:^:${n}\t:g" $i -i
done

# Combine all species data into single file
cat kra_tax/*_species.tsv > All_species.tsv

# Remove header duplicates
sed "/Species\t/d" All_species.tsv  -i 

###########################################################
# Step 11: Taxonomic Profiling with mOTUs
###########################################################
# Load mOTUs environment
source /datanode02/zhangzh/.apps/motu-env.sh

# Run mOTUs for metagenomic profiling (more accurate for microbes)
# -f/-r: Forward/reverse reads
# -q: Quiet mode
# -n: Sample name
# -c: Output coverage information
# -t 30: Threads
# >: Output profile
for i in `cat sample.name`; do
    motus profile  -f 02Cleandata/Trimgalore_out/${i}_1.fq.gz  \
        -r 02Cleandata/Trimgalore_out/${i}_2.fq.gz  \
        -q -n ${i} -c  -t 30 > 02Cleandata/reads_analysis/motu/Raw_result/${i}_motus.tsv
done 

cd 02Cleandata/reads_analysis/motu/

# Process mOTUs output for merging
# Extract abundance values
for i in Raw_result/* ; do
    awk -F"\t" 'NR>2{print $3}' ${i} > process/${i#*/}
done

# Extract taxonomic names and IDs
\ls Raw_result/*tsv | head -n 1 | xargs awk -F"\t" 'NR>2{print $1}'  > NAME
\ls Raw_result/*tsv | head -n 1 | xargs awk -F"\t" 'NR>2{print $2}'  > Tax

# Merge all samples into single abundance table
paste NAME process/* Tax > Abundance_merge.tsv

###########################################################
# Step 12: Taxonomic Profiling with SingleM
###########################################################
# Load SingleM environment
source /datanode02/zhangzh/.apps/singlem.sh

# Create symlink to SingleM metapackage (database)
ln -s /datanode02/zhangzh/database/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb 

# Run SingleM pipeline (16S rRNA gene based profiling)
# --metapackage: SingleM database package
# --threads 20: Threads
# -1/-2: Input paired-end reads
# -p: Output taxonomic profile
# --otu-table: Output OTU table
for i in `cat sample.name`; do
    singlem pipe --metapackage  S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb  \
        --threads 20    -1  01Rawdata/${i}_1.fq.gz    \
        -2 01Rawdata/${i}_2.fq.gz   -p  02Cleandata/reads_analysis/singlem/${i}_singlem.tsv   \
        --otu-table   02Cleandata/reads_analysis/singlem/${i}_otutab.tsv
done 

# Clean up symlink
rm ./S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb  -r 

# Combine all SingleM profiles into one file (keep header only once)
cat 02Cleandata/reads_analysis/singlem/*_singlem.tsv | sed 1!{/coverage/d}  > 02Cleandata/reads_analysis/Allsample.single.tsv 

cd 02Cleandata/reads_analysis

# Generate taxonomic summaries at different levels
# Phylum level
singlem summarise --input-taxonomic-profile Allsample.single.tsv \
    --output-species-by-site-relative-abundance Singlem.phylum_by_sample.tsv \
    --output-species-by-site-level phylum

# Class level
singlem summarise --input-taxonomic-profile Allsample.single.tsv \
    --output-species-by-site-relative-abundance Singlem.class_by_sample.tsv \
    --output-species-by-site-level class

# Genus level
singlem summarise --input-taxonomic-profile Allsample.single.tsv \
    --output-species-by-site-relative-abundance Singlem.genus_by_sample.tsv \
    --output-species-by-site-level genus

###########################################################
# Step 13: Clean Up Intermediate Files
###########################################################
# Move BAM files to central directory
for i in `cat sample.name`; do
    mv  05MAG/binning/${i}_matb2_maxb2_conc/work_files/${i}.bam  05MAG/bam
done 

# Remove large work directories
for i in `cat sample.name`; do
    rm -r 05MAG/binning/${i}_matb2_maxb2_conc/work_files
done 

###########################################################
# Step 14: Gene Expression Quantification (featureCounts/TPM)
###########################################################
# Load palmid environment (contains featureCounts dependencies)
source  /datanode02/zhangzh/.apps/palmid.sh

# Quantify gene expression with featureCounts
# -a: GFF annotation file
# -o: Output count file
# -t CDS: Feature type to count
# -g ID: Attribute to use as gene ID
# -T 30: Threads
# -p: Paired-end reads
# Last argument: BAM file with alignments
for i in `cat sample.name` ; do
    /datanode02/zhangzh/minisoft/subread-2.0.4-Linux-x86_64/bin/featureCounts \
        -a 04ORF/prodigal/${i}_contigs1k_prodigal.gff  \
        -o  04ORF/Counts/${i}_gene.counts \
        -t CDS -g ID -T 30  -p 05MAG/bam/${i}.bam
done 

# Remove BAM files to save space
rm 05MAG/bam/*.bam

# Calculate TPM (Transcripts Per Million) from raw counts
cd 04ORF
for  i in Counts/*counts; do
    n=$(basename $i .counts) ;
    /datanode02/zhangzh/minisoft/TPM_featureCounts.sh ${i} TPM/${n}.txt
done
cd ../

###########################################################
# Step 15: Functional Annotation with KOFAMscan (KEGG Orthologs)
###########################################################
# Load KOFAMscan environment
source /datanode02/zhangzh/.apps/kofam.sh
source /datanode02/zhangzh/.apps/DAS_tool.sh

# Run KOFAMscan for KEGG ortholog annotation
# -E 1e-5: E-value threshold
# --cpu 30: Threads
# -f detail-tsv: Output format (detailed TSV)
# -o: Output file
# Last argument: Protein sequences (faa)
for i in `cat sample.name`;  do
    /datanode02/zhangzh/minisoft/kofamscan -E 1e-5 --cpu 30  \
        -f detail-tsv -o 04ORF/Annotation/KEGG/Raw_result/${i}_kegg_raw.txt   \
        04ORF/prodigal/${i}_contigs1k_prodigal.faa
done 

cd 04ORF/Annotation/KEGG
mkdir process -p
mkdir merge_tpm -p 

# Process KOFAMscan output:
# 1. Filter by E-value (<1e-5)
# 2. Add sample name prefix
# 3. Sort by gene ID and E-value
# 4. Keep only best hit per gene
for i in  Raw_result/*_kegg_raw.txt; do
    n=$(basename $i _kegg_raw.txt);
    awk -F"\t" 'NR>2 && $6<1e-5 {print $0}' OFS="\t"  ${i} |
        sed "s:^:${n}\t:g" |
        sort -t$'\t' -s -k3n,3 -k7g,7r |
        awk -F"\t" '!a[$3]++'  OFS="\t"   > process/${n}.tsv
done

# Merge KEGG annotations with TPM values
for i in  Raw_result/*_kegg_raw.txt; do
    n=$(basename $i _kegg_raw.txt);
    awk -F"\t" 'NR==FNR{a[$3]=$0}NR>FNR&&a[b=$1]{print $0, a[b]}' OFS="\t"  \
        process/${n}.tsv  ../../TPM/${n}_gene.txt |
        awk -F"\t" 'BEGIN{print "ORFs\tCounts\tTPM\tSample\tKOs"}{print $1,$2,$4,$5,$8}' OFS="\t" > merge_tpm/${n}.tpm.tsv
done 

# Combine all samples into single KEGG abundance table
cat merge_tpm/*.tsv | sed '1!{/Sample/d}' > Allsample_KEGG_abundance.tsv 

# Reshape KEGG abundance table with R script (wide format)
source /datanode02/zhangzh/.apps/palmid.sh
Rscript /datanode02/zhangzh/Rfile/dcast-kegg.R Allsample_KEGG_abundance.tsv

###########################################################
# Step 16: MAG Coverage Calculation with CoverM
###########################################################
# Load CoverM environment
source /datanode02/zhangzh/.apps/coverm.sh

# Create temporary directory for CoverM
mkdir  coverm_tmp
export TMPDIR="coverm_tmp" 

# Calculate MAG coverage (TPM) with CoverM
# --coupled: Paired-end reads
# -t 30: Threads
# --methods tpm: Calculate TPM
# --genome-fasta-files: Input MAG fasta files
# -x fa: File extension for MAGs
# -o: Output file
for i in `cat sample.name`; do
    coverm genome --coupled 02Cleandata/Trimgalore_out/${i}_1.fq.gz  02Cleandata/Trimgalore_out/${i}_2.fq.gz \
        -t 30 --methods tpm \
        --genome-fasta-files 05MAG/MAG_data/MAG_95/dereplicated_genomes/* \
        -x fa -o 05MAG/MAG_data/Coverm/Raw_result/${i}_coverm.tsv
done 

# Clean up temporary directory
rm coverm_tmp -r 

cd 05MAG/MAG_data/Coverm/

# Process CoverM output for merging
for i in Raw_result/* ; do
    awk -F"\t" '{print $2}' ${i} > process/${i#*/}
done

# Extract MAG names
\ls Raw_result/*tsv | head -n 1 | xargs awk -F"\t" '{print $1}' OFS="\t" > NAME

# Merge all samples into single coverage table
paste NAME process/* > Abundance_merge.tsv

###########################################################
# Step 17: MAG Functional Annotation with KOFAMscan
###########################################################
# Load KOFAMscan environment
source /datanode02/zhangzh/.apps/kofam.sh
source /datanode02/zhangzh/.apps/DAS_tool.sh

# Move Prodigal output for MAGs to accessible directory
if ls "05MAG/MAG_data/prodigal"* &> /dev/null; then
    echo "prodigal file exist"
else
    echo "mv prodigal file to 05MAG/MAG_data"
    mv 05MAG/MAG_data/MAG_95/data/prodigal/  05MAG/MAG_data
fi

# Annotate MAG proteins with KOFAMscan
for i in 05MAG/MAG_data/MAG_95/dereplicated_genomes/* ; do
    /datanode02/zhangzh/minisoft/kofamscan -E 1e-5 --cpu 30  \
        -f detail-tsv -o 05MAG/MAG_data/Annotation/KEGG/Raw_result/${i##*/}_kegg_raw.txt   \
        05MAG/MAG_data/prodigal/${i##*/}.faa
done 

# Process MAG KEGG annotations:
# 1. Filter by E-value (<1e-5)
# 2. Add MAG name
# 3. Sort by gene ID and E-value
# 4. Keep only best hit per gene
# 5. Combine all MAG annotations
for i in 05MAG/MAG_data/Annotation/KEGG/Raw_result/* ; do
    n=$(basename $i _kegg_raw.txt) ;
    awk -F"\t" -v n="$n"   'NR>2 && $6<1e-5{print $2,$3,$5,$6,$7,n}' OFS="\t" ${i} |
        sort -t$'\t' -s -k1,1 -k4g,4r |
        awk -F"\t" '!a[$1]++'  OFS="\t"  >> 05MAG/MAG_data/Annotation/KEGG/MAG_kegg_uniq.txt
done

# Reshape MAG KEGG table with R script
source /datanode02/zhangzh/.apps/palmid.sh
cd 05MAG/MAG_data/Annotation/KEGG/
Rscript /datanode02/zhangzh/Rfile/dcast-kegg-mag.R MAG_kegg_uniq.txt 

###########################################################
# Step 18: MAG Quality Filtering and Taxonomic Classification
###########################################################
# Move Prodigal output to final location
mv 05MAG/MAG_data/MAG_95/data/prodigal/  05MAG/MAG_data 

# Extract high/medium quality MAG IDs (≥50% complete, ≤10% contaminated)
awk -F","  '$2>=50 && $3<=10{print $1}' 05MAG/MAG_data/MAG_95/data_tables/genomeInfo.csv > 05MAG/MAG_data/MHMAG.id

cd 05MAG/MAG_data

# Create batch file for GTDB-Tk (absolute paths)
sed "s:^:`pwd`\/DAS_bins/:g"  MHMAG.id  | awk -F"/" '{print $0,$NF}' OFS="\t" > gtdbtk_batch.tsv

# Load GTDB-Tk environment
source /apps/source/anaconda3.sh
source  activate /datanode02/zhangzh/.conda/envs/gtdbtk

# Classify MAGs with GTDB-Tk (species-level)
# --genome_dir: Directory with MAG fasta files
# --out_dir: Output directory
# --extension fa: File extension for MAGs
# --cpus 30: Threads
gtdbtk classify_wf --genome_dir MAG_95/dereplicated_genomes/  \
    --out_dir classify_wf_SGB  --extension fa  --cpus 30 > gtdb-sgb.log

# Infer phylogenetic trees for Archaea and Bacteria
gtdbtk infer --prefix arc  --gamma  \
    --msa_file classify_wf_SGB/align/gtdbtk.ar53.user_msa.fasta   \
    --out_dir classify_wf_SGB/infer --cpu 30

gtdbtk infer --prefix bac  --gamma  \
    --msa_file classify_wf_SGB/align/gtdbtk.bac120.user_msa.fasta   \
    --out_dir classify_wf_SGB/infer --cpu 30

###########################################################
# Step 19: MAG Functional Annotation with MetaCerberus
###########################################################
# Load MetaCerberus environment
source /datanode02/zhangzh/.apps/metacerberus.sh

# Annotate MAG proteins with multiple databases:
# KOFam_all (KEGG), COG, VOG, PHROG (phage), CAZy (carbohydrate-active enzymes)
# --protein: Protein sequences (faa)
# --hmm: Databases to use
# --dir_out: Output directory
# --cpus 30: Threads
for i in MAG_95/dereplicated_genomes/* ; do
    metacerberus.py  --protein  prodigal/${i##*/}.faa   \
        --hmm "KOFam_all, COG, VOG, PHROG, CAZy" \
        --dir_out Annotation/metacerberus/${i##*/} --cpus 30
done

# Optional: Alternative annotation for high/medium quality MAGs only
# for i in `cat  MHMAG.id`; do
#     metacerberus.py  --protein  prodigal/${i}.faa   \
#         --hmm "KOFam_all, COG, VOG, PHROG, CAZy" \
#         --dir_out Annotation/metacerberus/${i} --cpus 30
# done 
###########################################################
# Step 19: MAG Metabolic Potential Analysis with METABOLIC-G
###########################################################
# Load METABOLIC environment (conda)
source activate /datanode02/zhangzh/.conda/envs/metabolic4.sh

# Run METABOLIC-G for metabolic pathway prediction
# -in: Input directory with MAG protein sequences (Prodigal faa)
# -o: Output directory for metabolic analysis results
perl /datanode02/yut/Software/Metabolic/METABOLIC-G.pl \
    -in 05MAG/MAG_data/prodigal/ \
    -o 05MAG/MAG_data/Metabolic_re

echo "Pipeline completed successfully!"