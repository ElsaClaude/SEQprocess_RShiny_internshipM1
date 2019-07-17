##################################
######## Required lIBRARY ########
##################################
library(parallel)
library(GenomicRanges)
library(limma)
library(R.utils)
library(sequenza)
library(Biobase)
library(SummarizedExperiment)
# SEQprocess Report
library(data.table)
library(fastqcr)
library(pander)
library(knitr)
library(png)
library(grid)
library(gridExtra)
library(ggplot2)
library(reshape2)

##################################
### 1. configure paths & names ###
##################################
program.dir="/home/elsa/soft/miniconda3/envs/SEQprocess_internship"
program2.dir="/usr/bin"
reference.dir="/home/Tera/masterBI/stageM1/data"

# .idx
fq1.idx=".1.fastq$|_R1.fastq$|.1_val_1.fq$|_1.fq$|.1_val_1.fq$|.R1_val_1.fq$|R1_paired.fq$"
fq2.idx=".2.fastq$|_R2.fastq$|.2_val_2.fq$|_2.fq$|.2_val_2.fq$|.R2_val_2.fq$|R2_paired.fq$"
bam.idx=".rg.od.bam$|.rmdu.bam$|.realign.bam$|.recal.bam$|.Aligned.sortedByCoord.out.bam$|.sam$|accepted_hits.bam$"
vcf.idx=".f.vcf$|.muse_variants.vcf$|.snp.vcf$|.indel.vcf$|_variants.vcf$"
csv.idx=".hg38_multianno.csv$|.hg19_multianno.csv$"
# if using GDC pipeline variant call method, set up normal and tumor samples
n.sample="ASN"
t.sample="ASC"

##################################
######### 1.QualtiyCheck #########
##################################
# FASTQC
fastqc.path=file.path(program.dir, "bin/fastqc")

##################################
########### 2.Trimming ###########
##################################
# Trim_galore
trim_galore.path=file.path(program.dir, "bin/trim_galore")
trim.quality=30
trim.clip_R1=13
trim.clip_R2=13

# Cutadapt
cutadapt.path=file.path(program2.dir, "cutadapt")
m=17 # minimum read length
adapt.seq="TGGAATTCTCGGGTGCCAAGG"

##################################
########## 3. Alignment ##########
##################################
# BWA-mem & BWA-aln
bwa.path=file.path(program2.dir, "bwa")

# TopHat2
tophat2.path=file.path(program.dir, "bin/tophat2")

# Bowtie2
Sys.setenv("PATH"=paste(Sys.getenv("PATH"), "/home/elsa/soft", sep=":"))
bowtie2.path=file.path(program.dir, "bin/bowtie2")

# STAR
STAR.path=file.path(program.dir, "bin/STAR")
sjdbOverhang=100
outFilterMultimapScoreRange=1 
outFilterMultimapNmax=20
outFilterMismatchNmax=10
alignIntronMax=500000
alignMatesGapMax=1000000
sjdbScore=2
alignSJDBoverhangMin=1
outFilterMatchNminOverLread=0.33 
outFilterScoreMinOverLread=0.33 

# samtools 
samtools.path=file.path(program.dir, "bin/samtools")

# Picard AddOrReplaceReadGroups & Picard ReorderSam
picard.path=file.path(program.dir, "share/picard-2.20.2-0/picard.jar")
RGLB="LC"
RGPL="Illumina"
RGPU="runbarcode"
SORT_ORDER="coordinate"
VALIDATION_STRINGENCY="LENIENT"

ALLOW_INCOMPLETE_DICT_CONCORDANCE=FALSE 
ALLOW_CONTIG_LENGTH_DISCORDANCE=FALSE
CREATE_INDEX=TRUE

##################################
###### 4. Remove Duplicates ######
##################################
# Picard MarkDupliecates
REMOVE_DUPLICATES=TRUE
#VALIDATION_STRINGENCY="LENIENT"


##################################
######## 5. Re-Alignment #########
##################################
# GATK
GATK.path=file.path("/home/elsa/soft/GenomeAnalysisTK-3.7-0-gcfedb67/GenomeAnalysisTK.jar") # variant calling too

##################################
###### 6. Variant Calling ########
################################## # MuTect2 path=GATK path
# HaplotypeCaller
genotyping_mode="DISCOVERY"
output_mode="EMIT_VARIANTS_ONLY"
stand_call_conf_number=30
FS=30.0
QD=2.0
FS=trimws(format(round(FS, 1), nsmall = 1))
QD=trimws(format(round(QD, 1), nsmall = 1))
QUAL=50
DP=5
gatk.window=35
cluster=3

# VarScan2
varscan.path=file.path(program.dir, "share/varscan-2.4.3-2/VarScan.jar")
mapQ=1
min_tumor_freq=0.1
max_normal_freq=0.05
p_value=0.07

# MuSE
MuSE.path=file.path(program.dir, "bin/MuSE")
#bgzip.path
#tabix.path

# MuTect2
minN=2
filteredrecordsmergetype="KEEP_IF_ANY_UNFILTERED"
contamination_fraction_to_filter=0.02

# SomaticSniper
somaticsniper.path=file.path(program.dir, "bin/bam-somaticsniper")
mapQual=1
LOH=TRUE
Genotype=TRUE
somaticQual=15
somaticMutation=0.01
Theta=0.85
Hap.number=2
Hap.diff=0.001
out.format="vcf"

##################################
######## 7. Annotation ###########
##################################
# VEP
vep.path=file.path(program.dir, "bin/vep_install.pl")
vep.db.dir=file.path(reference.dir, "ngs_ref/vep")
perl5.10.path="/usr/bin"

# ANNOVAR
annovar.db.dir=file.path("home/elsa/soft/annovar/humandb")
vcf2annovar.pl=file.path("home/elsa/soft/annovar/convert2annovar.pl")
format = "vcf4"
coverage = 0

table_annovar.pl=file.path("home/elsa/soft/annovar/table_annovar.pl")
ref = "hg38"
protocol = "knownGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,exac03,avsnp147,ljb26_all,cosmic70"
protocol.type = "g,r,r,f,f,f,f,f,f"
nastring = "."

##################################
####### 8. RNA-abundance #########
##################################
# Cufflinks
cufflinks.path=file.path(program.dir, "bin/cufflinks")

# HT-Seq
htseq.path=file.path(program.dir, "bin/htseq-count")
Mode="intersection-nonempty" # read count method
stranded="no" #For stranded=no, a read is considered overlapping with a feature regardless of whether it is mapped to the same or the opposite strand as the feature.
idattr="gene_id" #GFF attribute to be used as feature ID.
htseq.r="pos" #For pos, read alignments whose mate alignment have not yet been seen are kept in a buffer in memory until the mate is found.
htseq.a=10 #skip all reads with alignment quality lower than the given minimum value


##################################
######## 9. Copy number ##########
##################################
# Sequenza
sequenza.util=system.file("exec", "sequenza-utils.py", package="sequenza")
window=1000000


##################################
######### 10. Make Set ###########
##################################
# GATK DepthOfCoverage
minBaseQuality=1
minMappingQuality=1


#####################################################
# Reference files  for alignmnet and variant calling 
#####################################################
#genome fasta file #Version = GRCH38
ref.fa=file.path(reference.dir,"GRCh38_no_alt_analysis_set.201503031.fa")
chrom.fa= file.path(reference.dir, "chroms/")

## Reference index #######
# bwa
bwa.idx=file.path(reference.dir, "GRCh38_no_alt_analysis_set.201503031")
# bowtie2
bowtie.idx=file.path(reference.dir, "GRCh38_no_alt_analysis_set.201503031")
# STAR
#star.idx.dir=file.path(reference.dir, "ngs_ref/star/gencode.v22")
star.idx.dir=file.path(reference.dir,"ngs_ref/star/gencode.v27")
# tophat2
transcriptome.idx=file.path(reference.dir, "ngs_ref/ftp.ebi.ac.uk/pub/databases/blueprint/reference/20150407_reference_files/transcript_index/GRCH38")

# Reference files for variant Calling
ref.gold_indels=file.path(reference.dir, "resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf")
ref.dbSNP=file.path(reference.dir, "common_all_20180418.vcf")
fn.hapmap.vcf=file.path(reference.dir, "resources_broad_hg38_v0_hapmap_3.3.hg38.vcf")
fn.omni.vcf=file.path(reference.dir, "resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf")
fn.1000g.vcf=file.path(reference.dir, "resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf")
#cosmic.vcf=file.path(reference.dir, "COSMIC/CosmicCodingMuts_v76.vcf") 

#####################################################
# Reference files  for RNAquantitiation
#####################################################
# gtf, gff 
#ref.gtf=file.path(reference.dir, "ngs_ref/ftp.ebi.ac.uk/pub/databases/blueprint/reference/20150407_reference_files/gencode.v22.annotation.201503031.gtf")
ref.gtf=file.path(reference.dir, "gencode.v30.annotation.gtf")
mir.gff=file.path(reference.dir, "gencode.v30.chr_patch_hapl_scaff.annotation.gff3")
refGene.path=file.path(reference.dir, "ucsc_table/refGene_hg38")