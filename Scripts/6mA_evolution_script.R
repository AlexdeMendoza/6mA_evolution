library(data.table)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(bsseq)
library(R.utils)
library(reshape2)
library(ggseqlogo)
library(cowplot)

################################################################################
## Load genomic files and get 6mA in genes and in nucleotides
################################################################################

# Name the sample
name <- "Cfra_Multinucleated"

# Establish the order of species for the plots
All_species_order <- c("Cfra", "Aapp", "Awhi", "Cper", "Clim", "Spun", "Acas", "Ddis", "Crei", "Cvar", "Otau", "Mpus", "Cpar", "Ptri", "Alim", "Ngru", "Nfow", "Tvag")
Species_order <- c("Cfra", "Aapp", "Awhi", "Cper", "Clim", "Spun", "Acas", "Crei", "Cvar", "Otau", "Mpus", "Alim", "Ngru", "Nfow", "Tvag")


# ------------------------------------------------------------------------------
# Load genomic information and filtrate file
# ------------------------------------------------------------------------------


#########Unix command line code to produce a genomic file from nanopore#########
#
## Run nanopore_guppy to do the basecalling 
#
#nanopore_guppy guppy_basecaller -i fast5/ -s guppy_GPU_allmods_R9/ -d /rerio/basecall_models/ -c res_dna_r941_min_modbases-all-context_v001.cfg --device auto --bam_out --recursive --align_ref Genome.fasta --compress_fastq
#
## Run samtools to Merge BAM files into a single BAM file
#samtools merge -@ 20 Allmods.merged.bam *bam
#
## Run samtools to sort the merged BAM file
#samtools sort -@ 20 Allmods.merged.bam -o Allmods.merged.sorted.bam
#
## Run samtools to index the sorted BAM file
#samtools index -@ 20 Allmods.merged.sorted.bam
#
## Run modbam2bed to convert the sorted BAM file to BED format
#modbam2bed -e -m 6mA -t 30 Genome.fasta guppy_GPU_allmods_R9/Allmods.merged.sorted.bam > Species.R9_6mA.methyl.bed 
#modbam2bed -e -m 5mC -t 30 Genome.fasta guppy_GPU_allmods_R9/Allmods.merged.sorted.bam > Species.R9_5mC.methyl.bed
#
## Run bedtools to create coordinates with neighbouring sequence context for every position
# bedtools slop -l 1 -r 2 -s -i Species.R9_6mA.methyl.bed -g Genome.fasta.fai > Species.R9_6mA.methyl.plus3.bed
#
## Run bedtools to extract sequence context from the genome FASTA file based on the coordinates 
# bedtools getfasta -fi Genome.fasta -fo Species.R9_6mA.methyl.plus3.fa -s -bed Species.R9_6mA.methyl.plus3.bed
#
## Run seqkit to convert FASTA format to tab-separated format and extracting sequences
# seqkit fx2tab Species.R9_6mA.methyl.plus3.fa | cut -f 2 > a
#
# Merge BED file containing genomic coordinates with the extracted sequences
# paste Species.R9_6mA.methyl.bed a > Species.R9_6mA.methyl.context.bed
#
################################################################################

# Cfrag.R9_6mA.methyl.context.bed.gz is used as example, but all methylBeds were treated the same way.

# Load the .bed file with the genomic information generated before
df <- fread("../Cfrag.R9_6mA.methyl.context.bed.gz")

# Get the % of each context
AT <- df %>% filter(grepl(V17, pattern = "[ATCG]AT[ATCG]"))
AA <- df %>% filter(grepl(V17, pattern = "[ATCG]AA[ATCG]"))
AC <- df %>% filter(grepl(V17, pattern = "[ATCG]AC[ATCG]"))
AG <- df %>% filter(grepl(V17, pattern = "[ATCG]AG[ATCG]"))

# Calculate the global rates for the 4 combinations
global_stats <- data.frame(mAT =100*sum(AT$V13)/(sum(AT$V12)+sum(AT$V13)),
                           mAA =100*sum(AA$V13)/(sum(AA$V12)+sum(AA$V13)),
                           mAC =100*sum(AC$V13)/(sum(AC$V12)+sum(AC$V13)),
                           mAG =100*sum(AG$V13)/(sum(AG$V12)+sum(AG$V13)))

rm(AA,AC,AG)

# Extract AT information from genomic data
dat <- AT %>% dplyr::rename(chr = V1, start = V2, end = V3, strand = V6, mC_perc = V11,
                            C_reads = V13) %>% mutate(CT_reads = V12 + C_reads)

# Generate a BS object, pretending that ApTs are CpGs
bs_obj <- BSseq(gr = GRanges(seqnames = dat$chr,
                             ranges = IRanges(start = as.numeric(dat$start+1),
                                              end = as.numeric(dat$start+1)), strand = dat$strand ), 
                sampleNames = name, M = as.matrix(dat$C_reads), Cov = as.matrix(dat$CT_reads),
                rmZeroCov = FALSE)

# Get all methylation levels per all ApT sites and see what is the distribution of values
getMeth(bs_obj, type = "raw", what = "perBase") %>% hist(breaks = 100, main = "Distribution of 6mA levels", xlab = "6mA %")



# ------------------------------------------------------------------------------
# Plot the ApT strand correlation values of all species
# ------------------------------------------------------------------------------

# Calculate correlation between the methylation in + and - strands 

AT_dinuc <- dat %>% mutate(locus = if_else(strand == "+", paste0(chr,":",start+1),paste0(chr,":",start) )) %>%
  dplyr::select(locus, strand, mC_perc, CT_reads)

# Filter positions with 10x coverage in both strands

AT_dinuc_watson <- AT_dinuc %>% filter(CT_reads >= 10, strand == "+")
AT_dinuc_crick <- AT_dinuc %>% filter(CT_reads >= 10, strand == "-")

# Take only those present in both data frames (10x coverage in both strands)
AT_dinuc <- inner_join(AT_dinuc_watson, AT_dinuc_crick, by = "locus")

# Visualize and compute the correlation
smoothScatter(AT_dinuc$mC_perc.x,AT_dinuc$mC_perc.y)
cor(AT_dinuc$mC_perc.x,AT_dinuc$mC_perc.y, method = "pearson")
cor(AT_dinuc$mC_perc.x,AT_dinuc$mC_perc.y, method = "spearman")

# Merge the values from + and - strand into a single "ApT" value
bs_obj_collapsed <- strandCollapse(bs_obj)

# Save the BSobject as a RDS object,
saveRDS(object = bs_obj_collapsed, file = paste0(name,".AT_bsseq.rds"))

# Calculate Pearson and Spearman correlation

AT_pearson <- cor(AT_dinuc$mC_perc.x,AT_dinuc$mC_perc.y, method = "pearson")

AT_spearman <- cor(AT_dinuc$mC_perc.x,AT_dinuc$mC_perc.y, method = "spearman")

# Compute the correlation values for all species
Cfra_ATCorrelation <- c(name, AT_pearson, AT_spearman )#Repeat with the rest of species


# Prepare a matrix with all species values
AT_Pearson_Correlation <- matrix(
  c(
    Cfra_ATCorrelation[1], Cfra_ATCorrelation[2],
    Aapp_ATCorrelation[1], Aapp_ATCorrelation[2],
    # ... Repeat for the rest of species
    Tvag_ATCorrelation[1], Tvag_ATCorrelation[2]
  ),
  ncol = 2,
  byrow = TRUE
)

# Prepare the dataframe for plot generation
AT_Pearson_Correlation_DF <- as.data.frame(AT_Pearson_Correlation)
colnames(AT_Pearson_Correlation_DF) <- c("Species", "Pearson")
AT_Pearson_Correlation_DF_long <- tidyr::gather(AT_Pearson_Correlation_DF, key = "variable", value = "value", -Species)
AT_Pearson_Correlation_DF_long$value <- as.numeric(AT_Pearson_Correlation_DF_long$value)


# Generate the bar chart with correlation values of each species
AT_Pearson_Correlation_Bar_Chart<- ggplot(AT_Pearson_Correlation_DF_long, aes(x = reorder(Species, match(Species, Species_order)), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(title = "CpG symmetric methylation correlation", x = 'Species', y = "Pearson correlation value", fill = "Variable") +
  theme_light() +
  guides(fill = "none") +
  scale_y_continuous(limits = c(min(df_longA$value), max(df_longA$value)))



# ------------------------------------------------------------------------------
# Load gene bed file and get the 6mApT levels on gene bodies
# ------------------------------------------------------------------------------


# Function to import a 6 column bed file into R as a Genomic Ranges object

read_bed6_to_GRobject <- function(bedfile, chr_append = F){
  dat <- fread(bedfile)
  if (chr_append == T){
    dat$V1 <- paste0("chr",dat$V1)
  }
  #dat <- dat[dat$V1 %in% names(sLengths),]
  if (sum(grepl(unique(dat$V6), pattern = "C")) == 1){
    dat$V6 <- ifelse(dat$V6 == "C", yes = "-","+")
  }
  gr <- GRanges(seqnames = Rle(dat$V1),
                ranges = IRanges(start = dat$V2, end = dat$V3),
                strand = dat$V6, 
                gene_id = dat$V5, target_id = dat$V4)
  return(gr)
}

add_loci <- function(gr){
  
  loci <- paste(seqnames(gr), start(gr), sep=":") %>%
    paste(end(gr), sep = "-")
  
  gr$loci <- loci
  
  return(gr)
  
}

# all gene beds are "chr start end name ID strand" format.

# Load the gene bed file for a species
genes <- read_bed6_to_GRobject("PASA_Creolimax2.gene.bed") 

# Obtain mean coverage per sample
getCoverage(bs_obj_collapsed,type = "Cov") %>% colMeans()


# Function to extract methylation information for genomic ranges
get_mA_on_ranges <- function(bs_obj, gene_gr, thres = 0.1){
  coverage_positions <- getCoverage(bs_obj, regions = gene_gr, type = "Cov", what = "perBase")
  position_per_gene <- lapply(coverage_positions, FUN = length) %>% unlist()
  methylation_positions <- getMeth(bs_obj, regions = gene_gr, type = "raw", what = "perBase")
  meth_position_per_gene <- lapply(methylation_positions, function(x) sum(x > thres)) %>% unlist()
  
  mean_coverage_nanopore <- getCoverage(bs_obj, regions = gene_gr, type = "Cov", what = "perRegionAverage")[,1]
  M <- getCoverage(bs_obj, regions = gene_gr, type = "M", what = "perRegionTotal")
  C <- getCoverage(bs_obj, regions = gene_gr, type = "Cov", what = "perRegionTotal")
  mCG_on_nanopore <- (M/C)[,1]
  rm(M,C)
  
  df <- data.frame(gene = gene_gr$target_id,
                   mean_coverage = mean_coverage_nanopore,
                   ApTs = position_per_gene,
                   methylated_ApT = meth_position_per_gene,
                   loci = add_loci(gene_gr) %>% .$loci, 
                   mApT = mCG_on_nanopore)
  return(df)
}

# Calculate the 6mAT (regional and per ApT) for every gene
genes_m6A <- get_mA_on_ranges(bs_obj = bs_obj_collapsed, gene_gr = genes)


# Check the coverage of each gene 
genes_m6A$mean_coverage %>% hist( breaks = 50, main = "Mean coverage per gene")

# Discard genes with mean coverage < 10x (this value was adjusted for species with low coverage)
total_gene_num <- genes_m6A %>% nrow() 
total_gene_covered <- genes_m6A %>% filter(mean_coverage >= 10) %>% nrow() 

# Filter those genes with 10x
genes_m6A_filt <- genes_m6A %>% filter(mean_coverage >= 10)

# Check the distribution of gene body methylation
genes_m6A_filt$mApT %>% hist( breaks = 50, main = "6mApT over gene bodies", xlab = "6mA %")

# Check how many genes have at least 10% methylation 
ten_perc_methylated <- genes_m6A_filt %>% filter(mApT >= 0.1) %>% nrow()

# Check the  distribution of "presence" of ApT per gene body
genes_m6A_filt$methylated_ApT %>% hist( breaks = 100, main = "Methylated ApTs  over gene bodies", xlab = "number of methylated ApTs per gene")

# Check how many genes have at least 3 methylated ApTs
genes_with_3_ApT <- genes_m6A_filt %>% filter( methylated_ApT >= 3) %>% nrow()


# ------------------------------------------------------------------------------
# Calculate number and percentage of methylated genes 
# ------------------------------------------------------------------------------


# Reneame the variable
Methylated_genes <- genes_with_3_ApT #From now on we asume a methylated gene has 3 ApTs more than a 10% of the time and has a coverage of at least 10X

# Calculate the numeber of unmethylated genes
Unmethylated_genes <- c(total_gene_covered- Methylated_genes )

# Calculate the percentage of unmethylated genes
Percentage_of_unmethylated_genes <- c(100*(Unmethylated_genes/total_gene_covered))

# Calculate the percentage of methylated genes
Percentage_of_methylated_genes <- c(100*(Methylated_genes/total_gene_covered))


# Create dataframe with the with the amounts of selected genes, methylated genes and unmethylated genes
Cfra_Multinucleated_Methylated_vs_unmethylated_genes <- data.frame(name, genes, total_gene_covered, Methylated_genes, Unmethylated_genes, Percentage_of_unmethylated_genes, Percentage_of_methylated_genes) #Repeat for all species


# Create a dataframe with the amounts of selected genes, methylated genes and unmethylated genes of all species
Methylated_vs_unmethylated_genes <- data.frame(Cfra_Multinucleated_Methylated_vs_unmethylated_genes, Apar_Methylated_vs_unmethylated_genes, Awhi_Methylated_vs_unmethylated_genes, Cperk_Methylated_vs_unmethylated_genes, Clim_Methylated_vs_unmethylated_genes, Spun_Methylated_vs_unmethylated_genes, Acas_Methylated_vs_unmethylated_genes,
                                               Crei_Methylated_vs_unmethylated_genes, Cvar_Methylated_vs_unmethylated_genes, Otau_Methylated_vs_unmethylated_genes, Mpus_Methylated_vs_unmethylated_genes, Alim_Methylated_vs_unmethylated_genes, Ngru_Methylated_vs_unmethylated_genes, Nfow_Methylated_vs_unmethylated_genes, Tvag_Methylated_vs_unmethylated_genes )
# Transpose dataframe
Methylated_vs_unmethylated_genes <- transpose(Methylated_vs_unmethylated_genes)

#Prepare the dataframes for the plot generation
Percentage_of_genes_data_long <- gather(Methylated_vs_unmethylated_genes, key = "GeneStatus", value = "Percentage",
                               Percentage_of_unmethylated_genes, Percentage_of_methylated_genes)
Number_of_genes_data_long <- gather(Methylated_vs_unmethylated_genes, key = "GeneStatus", value = "NumberOfGenes",
                                    Unmethylated_genes, Methylated_genes)

# Plot the percentage of methylated genes
Percentage_Genes_Methylation_Plot <- ggplot(Percentage_of_genes_data_long, aes(x = reorder(Species, match(Species, Species_order)), y = Percentage, fill = reorder(GeneStatus, -Percentage))) +
  geom_bar(stat = "identity", position = 'stack') +
  geom_text(aes(label = scales::percent(Percentage, scale = 1)), position = position_stack(vjust = 0.5), color = "black", size = ) +
  labs(title = "Percentage of Methylated vs Unmethylated Genes",
       x = "Species",
       y = "Genes Percentage") +
  scale_fill_manual(values = c("#fc8d62", "#66c2a5"), name = "Gene Status",
                    breaks = c('Percentage_of_unmethylated_genes', 'Percentage_of_methylated_genes'),
                    labels = c("Unmethylated Genes", "Methylated Genes")) +
  theme_minimal()


# Plot the number of methylated genes
Number_Genes_Methylation_Plot <- ggplot(Number_of_genes_data_long, aes(x = reorder(Species, match(Species, Species_order)), y = NumberOfGenes, fill = reorder(GeneStatus, -NumberOfGenes))) +
  geom_bar(stat = "identity", position = 'stack') + 
  geom_text(aes(label = NumberOfGenes), position = position_stack(vjust = 0.5), color = "black", size = 3) +
  labs(title = "Number of Methylated vs Unmethylated Genes",
       x = "Species",
       y = "Number of Genes") +
  scale_fill_manual(values = c("#fc8d62", "#66c2a5"), name = "Gene Status",
                    breaks = c('Unmethylated_genes', 'Methylated_genes'),
                    labels = c("Unmethylated Genes", "Methylated Genes")) +
  theme_minimal()




# Save methylated genes
write.table(genes_m6A, file = "Cfrag_Multinucleated_genes_6mA.tsv",
            quote = FALSE, sep = "\t",row.names = F)



# ------------------------------------------------------------------------------
# Calculate 6mA percentages in different nucleotide context
# ------------------------------------------------------------------------------


######Unix command line code to obtain methylation levels in different nucleotide contexts from a bed file######
#
## Calculate total coverage, methylated counts and methylation percentage in different nucleotide context from a bed file with genomic methylation information
## Print the total coverage, total methylated count, and percentage of methylation to a file 
#
#cat Species.R9_6mA.methyl.bed | awk '{cov+=$12; mA+=$13;} END{print "6mA", cov+mA, mA, mA/(cov+mA)*100;}' > a
#cat Species.R9_5mC.methyl.context.bed | awk '$17 ~ "[ATCG]CG[ATCG]"' | awk '{cov+=$12; mCG+=$13;} END{print "mCG", cov+mCG, mCG, mCG/(cov+mCG)*100;}' > b
#cat Species.R9_5mC.methyl.bed | awk '{cov+=$12; mC+=$13;} END{print "5mC", cov+mC, mC, mC/(cov+mC)*100;}' > c
#cat Species.R9_6mA.methyl.context.bed | awk '$17 ~ "[ATCG]AT[ATCG]"' | awk '{cov+=$12; mAT+=$13;} END{print "mAT", cov+mAT, mAT, mAT/(cov+mAT)*100;}' > d
#cat Species.R9_6mA.methyl.context.bed | awk '$17 ~ "[ATCG]AA[ATCG]"' | awk '{cov+=$12; mAA+=$13;} END{print "mAA", cov+mAA, mAA, mAA/(cov+mAA)*100;}' > e
#cat Species.R9_6mA.methyl.context.bed | awk '$17 ~ "[ATCG]AC[ATCG]"' | awk '{cov+=$12; mAC+=$13;} END{print "mAC", cov+mAC, mAC, mAC/(cov+mAC)*100;}' > f
#cat Species.R9_6mA.methyl.context.bed | awk '$17 ~ "[ATCG]AG[ATCG]"' | awk '{cov+=$12; mAG+=$13;} END{print "mAG", cov+mAG, mAG, mAG/(cov+mAG)*100;}' > g
#cat Species.R9_6mA.methyl.context.bed | awk '$17 ~ "[CG]AT[CG]"' | awk '{cov+=$12; SmATS+=$13;} END{print "SmATS", cov+SmATS, SmATS, SmATS/(cov+SmATS)*100;}' > h
#cat Species.R9_6mA.methyl.context.bed | awk '$17 ~ "[AT]AT[AT]"' | awk '{cov+=$12; WmATW+=$13;} END{print "WmATW", cov+WmATW, WmATW, WmATW/(cov+WmATW)*100;}' > i
#cat Species.R9_6mA.methyl.context.bed | awk '$17 ~ "[ACG]AT[CGT]"' | awk '{cov+=$12; VmATB+=$13;} END{print "VmATB", cov+VmATB, VmATB, VmATB/(cov+VmATB)*100;}' > j
#
# Create a header for the output file
#
#echo "context coverage m m_perc" > aa
#
# Concatenate all result files into one file
#
#cat aa a b c d e f g h i j > Species_R9_methyl_stats
#
################################################################################

# input the data from the previous bash script into R
df <- fread("SpeciesMethStats.txt")

df$Species <- factor(df$Species, levels = df$Species)

mALevels <- ggplot(data = df, aes(x = reorder(Species, match(Species, All_species_order)), y = mA)) +
  geom_bar(stat = "identity", fill="steelblue") +
  labs(y = "percentage 6mA (%)") +
  theme_light()

mAALevels <- ggplot(data = df, aes(x = reorder(Species, match(Species, All_species_order)), y = mAA)) +
  geom_bar(stat = "identity", fill='steelblue') +
  labs(y = "percentage mAA (%)") +
  theme_light()

mATLevels <- ggplot(data = df, aes(x = reorder(Species, match(Species, All_species_order)), y = mAT)) +
  geom_bar(stat = "identity", fill='steelblue') +
  labs(y = "percentage mAT (%)") +
  theme_light()

mACLevels <- ggplot(data = df, aes(x = reorder(Species, match(Species, All_species_order)), y = mAC)) +
  geom_bar(stat = "identity", fill='steelblue') +
  labs(y = "percentage mAC (%)") +
  theme_light()

mAGLevels <- ggplot(data = df, aes(x = reorder(Species, match(Species, All_species_order)), y = mAG)) +
  geom_bar(stat = "identity", fill='steelblue') +
  labs(y = "percentage mATG (%)") +
  theme_light()

SmATSLevels <- ggplot(data = df, aes(x = reorder(Species, match(Species, All_species_order)), y = SmATS)) +
  geom_bar(stat = "identity", fill='steelblue') +
  labs(y = "percentage SmATS (%)") +
  theme_light()

WmATWLevels <- ggplot(data = df, aes(x = reorder(Species, match(Species, All_species_order)), y = WmATW)) +
  geom_bar(stat = "identity", fill="steelblue") +
  labs(y = "percentage WmATW (%)") +
  theme_light()

VmATBLevels <- ggplot(data = df, aes(x = reorder(Species, match(Species, All_species_order)), y = VmATB)) +
  geom_bar(stat = "identity", fill="steelblue") +
  labs(y = "percentage VmATB (%)") +
  theme_light()

combo_plots <- plot_grid(mALevels, mATLevels, mAALevels, mACLevels, mAGLevels, SmATSLevels, WmATWLevels, VmATBLevels, ncol = 1 )



################################################################################
### Comparison between 6mAT methylation and transcription
################################################################################


#########Unix command line code to produce an abundance transcriptomic file from raw RNA-seq fastq file#########
#
## Building HISAT2 index using the genome FASTA file
#hisat2-build Genome.fasta Genome.hisat2_index
#
## Assigning FASTQ file names to variables
#fastq1=Species_RNA_seq_1.fastq
#fastq2= Species_RNA_seq_2.fastq
#
## Extracting file name without the ".fastq" extension
#name=$(echo $fastq1 | perl -pe "s/.fastq//")
#
## Run hisat2 to align paired-end reads
#hisat2 --max-intronlen 40000 --dta -N 1 -p $NSLOTS --summary-file ${name}.hisat2.log -x Genome.hisat2_index -1 $fastq1 -2 $fastq2 -S ${name}.hisat2.sam
#
##(For species with one fastq file) 
##(hisat2 --max-intronlen 40000 --dta -N 1 -p $NSLOTS --summary-file ${name}.hisat2.log -x Acas.hisat2_index -U $fastq1 -S ${name}.hisat2.sam)
#
##Run sambamba to convert SAM file to BAM format
#sambamba view -t $NSLOTS -f bam -S -o ${name}.hisat2.bam ${name}.hisat2.sam
#
##Run sambamba to sort the BAM file
#sambamba sort -t $NSLOTS -o ${name}.hisat2.sorted.bam ${name}.hisat2.bam
#
#
##Run stringtie to assemble transcripts (if the library is stranded)
#stringtie --rf -G Genomic.gtf -A Species.RNAseq_abundance.tsv -e ${name}.hisat2.sorted.bam
#
################################################################################


# Load RNA-seq data
tpm <-fread("Cfra_Multinucleated.RNAseq_abundance.tsv") %>% dplyr::rename(gene = "Gene\ ID") %>% 
  dplyr::select(gene, TPM)

# Join RNA and 6mA data
genes_m6A_filt <- inner_join(genes_m6A_filt, tpm)

# Check correlation between 6mA and TPM
cor(genes_m6A_filt$mApT, genes_m6A_filt$TPM, method = "spearman") 
cor(genes_m6A_filt$methylated_ApT, genes_m6A_filt$TPM, method = "spearman", use = "complete.obs")


# Assign decile category to genes based on their TPM values
genes_m6A_filt$decile_tpm <- cut(genes_m6A_filt$TPM, breaks = c(0,quantile(genes_m6A_filt$TPM[genes_m6A_filt$TPM >= 1], seq(0, 1, length = 11), type = 5)),
                                 labels = c("No expr.","0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))


# Adjust the decile category for genes with TPM = 0 to "No expr."
genes_m6A_filt <- genes_m6A_filt %>% mutate(decile_tpm = if_else(TPM == 0,"No expr.", decile_tpm )  )

# Convert the decile_tpm variable to a factor with specified levels
genes_m6A_filt$decile_tpm <- factor(genes_m6A_filt$decile_tpm, levels = c("No expr.","0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))

# Get unique decile values for each species
Cfra_Multinucleated_decile_values <- unique(genes_m6A_filt$decile_tpm)

# Loop through each unique decile value to extract gene lists for each decile
for (decile_value in Cfra_Multinucleated_decile_values) {
  # Subset data for the current decile value
  subset_data <- genes_m6A_filt[genes_m6A_filt$decile_tpm == decile_value, ]
  # Write the subset data to a file
  write.table(subset_data, file = paste0("Ngru_", gsub("\\s+", "", decile_value), ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}

# Assign a decile category to genes based on their 6mA levels
genes_m6A_filt$decile_6mA <- ntile(genes_m6A_filt$mApT, 10)


# Create a boxplot comparing gene body 6mApT % with gene expression deciles
gg_6mA_vs_tx_genes <- ggplot(genes_m6A_filt, aes(x = decile_tpm, y = 100*mApT)) + geom_boxplot(outlier.shape = NA) + 
  theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90),legend.position="none") + ylim(c(0,8)) +
  theme_classic() + ylab("gene body 6mApT %") + xlab("Gene expression decile (low to high)") +
  ggtitle(paste0("Creolimax gene body methylation vs transcription (n=",nrow(genes_m6A_filt),")"))

# Create a boxplot comparing the number of methylated mApT sites with gene expression deciles
gg_6mAsites_vs_tx_genes <- ggplot(genes_m6A_filt, aes(x = decile_tpm, y = methylated_ApT)) + geom_boxplot(outlier.shape = NA) + 
  theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90),legend.position="none") + ylim(c(0,10)) +
  theme_classic() + ylab("Methylated mApT sites") + xlab("Gene expression decile (low to high)") +
  ggtitle(paste0("Creolimax gene body methylation vs transcription (n=",nrow(genes_m6A_filt),")"))

# Create a boxplot comparing gene expression with 6mA deciles
gg_tx_genes_vs_6mA <- ggplot(genes_m6A_filt, aes(x = as.factor(decile_6mA), y = log10(TPM+1))) + geom_boxplot(outlier.shape = NA) + 
  theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90),legend.position="none") + 
  theme_classic() + ylab("log10(TPM+1)") + xlab("Gene body 6mA decile (low to high)") +
  ggtitle(paste0("Creolimax gene body methylation vs transcription (n=",nrow(genes_m6A_filt),")"))

# Combine individual plots into a single grid
combo_plots <- plot_grid(gg_6mA_vs_tx_genes,gg_6mAsites_vs_tx_genes, gg_tx_genes_vs_6mA, ncol = 1 )



################################################################################
### Gene orientation analysis and impact on 6mAT patterns
################################################################################

# Function to classify genes based on upstream gene direction
classify_genes_and_calc_intergenic <- function(granges){
  # Extract information from the GRanges object
  a <- data.frame(granges)[1:(length(granges)-1),]
  b <- data.frame(granges)[2:length(granges),]
  colnames(a) <- paste0(colnames(a),"_a")
  colnames(b) <- paste0(colnames(b),"_b")
  
  df <- cbind(a,b)
  
  # Separate genes on the negative and positive strands
  neg_strand <- df %>% filter(strand_a == "-")
  plus_strand <- df %>% filter(strand_b == "+")
  
  # Classify gene orientation as "Head_to_Tail," "Head_to_Head," or "Unk" (Unknown)
  neg_strand <- neg_strand %>% mutate(orientation = ifelse(seqnames_a == seqnames_b & strand_b == strand_a, "Head_to_Tail",
                                                           ifelse(seqnames_a == seqnames_b & strand_b != strand_a, "Head_to_Head","Unk")),
                                      distance = start_b - end_a) %>% 
    dplyr::select(target_id_a, orientation, distance) %>%
    dplyr::rename(gene = target_id_a)
  plus_strand <- plus_strand %>% mutate(orientation = ifelse(seqnames_a == seqnames_b & strand_b == strand_a, "Head_to_Tail",
                                                             ifelse(seqnames_a == seqnames_b & strand_b != strand_a, "Head_to_Head","Unk")),
                                        distance = start_b - end_a) %>% 
    dplyr::select(target_id_b, orientation, distance) %>%
    dplyr::rename(gene = target_id_b)
  
  df <- rbind(neg_strand, plus_strand)
  return(df)
}

# Classify genes based on upstream gene direction and calculate intergenic distances
gene_class <- classify_genes_and_calc_intergenic(granges = genes) 

# Save the gene classification results to a file
write.table(gene_class, file = "Cfra_genes_orientation.tsv",
            quote = FALSE, sep = "\t",row.names = F)


# Generate bed files with all the combinations (Methylated and orientation)
genes_to_subset  <- data.frame(genes) %>% dplyr::select(seqnames, start, end, target_id, strand) %>% 
  dplyr::mutate(gene = target_id) %>% inner_join(gene_class) %>% 
  dplyr::select(seqnames,start,end,gene,target_id,strand,orientation) %>% 
  inner_join(genes_m6A_filt) %>%
  mutate(status = ifelse(methylated_ApT >= 3, "Methylated","Unmethylated")) %>%
  dplyr::select(-mean_coverage, -ApTs, -methylated_ApT,-loci,-mApT,-TPM, -decile_tpm, -decile_6mA)


# Create subsets based on gene methylation status
Methylated_genes <- genes_to_subset %>% filter(status == "Methylated") %>% .[,(1:6)]
Unmethylated_genes <- genes_to_subset %>% filter(status == "Unmethylated") %>% .[,(1:6)]

# Save the bed files for further analysis
write.table(Methylated_genes, file = "Cfrag_Multinucleated.Methyl.bed",quote = FALSE, sep = "\t",row.names = F, col.names = F)
write.table(Unmethylated_genes, file = "Cfrag_Multinucleated.Unmethyl.bed",quote = FALSE, sep = "\t",row.names = F, col.names = F)

# Create subsets based on gene orientation and methylation status
h2h_m <- genes_to_subset %>% filter(orientation == "Head_to_Head", status == "Methylated") %>% .[,(1:6)]
h2h_u <- genes_to_subset %>% filter(orientation == "Head_to_Head", status == "Unmethylated") %>% .[,(1:6)]
h2t_m <- genes_to_subset %>% filter(orientation == "Head_to_Tail", status == "Methylated") %>% .[,(1:6)]
h2t_u <- genes_to_subset %>% filter(orientation == "Head_to_Tail", status == "Unmethylated") %>% .[,(1:6)]


# Save the bed files for further analysis
write.table(h2h_m, file = "Cfrag_Multinucleated.head_to_head.methyl.bed",quote = FALSE, sep = "\t",row.names = F, col.names = F)
write.table(h2h_u, file = "Cfrag_Multinucleated.head_to_head.unmethyl.bed",quote = FALSE, sep = "\t",row.names = F, col.names = F)
write.table(h2t_m, file = "Cfrag_Multinucleated.head_to_tail.methyl.bed",quote = FALSE, sep = "\t",row.names = F, col.names = F)
write.table(h2t_u, file = "Cfrag_Multinucleated.head_to_tail.unmethyl.bed",quote = FALSE, sep = "\t",row.names = F, col.names = F)

# these files were later plotted with DeepTools2 computeMatrix function
##



################################################################################
### Mean global methylation levels across genomic features
################################################################################

# Define promoter region in genes 
Promoters <- GenomicFeatures::promoters(genes, upstream = 500, downstream = 0)

# Calculate the 6mA (regional and per ApT) in the promoter region of every gene
Promoters_m6A <- get_mA_on_ranges(bs_obj, gene_gr = Promoters)

# Filtrate the promoters that have a coverage of less than 10x
Promoters_m6A_filt <- Promoters_m6A %>% filter(mean_coverage >= 10)

# Calculate the mean methylation of the promoter region
Methylation_Promoters <- mean(Promoters_m6A_filt[,6]) 


# Define Transcription Start Site (TSS) region in genes
TSS <- resize(genes, width = 500, fix = 'start')

# Calculate the 6mA (regional and per ApT) in the TSS of every gene
TSS_m6A <- get_mA_on_ranges(bs_obj, gene_gr = TSS)

# Filtrate the TSSs that have a coverage of less than 10x
TSS_m6A_filt <- TSS_m6A %>% filter(mean_coverage >= 10)

# Calculate the mean methylation of the TSS
Methylation_TSS <- mean(TSS_m6A_filt[,6])


# Function to obtain the Transcription End Site of genes
Get_TES <- function(gr, up = 0, down = 1, ...){
  
  # flip the strand of original GRanges
  rg <- GenomicRanges::invertStrand(gr)
  
  # extract TSS position of this flipped GRanges
  tesGr <- GenomicRanges::promoters(x = rg, upstream = up, downstream = down, ...)
  
  # reset the strand
  tesGr <- GenomicRanges::invertStrand(tesGr)
  
  return(tesGr)
}

# Define Transcription End Site region in genes
TES <- Get_TES(genes, up = 250, down = 250)

# Calculate the 6mA (regional and per ApT) in the TES of every gene
TES_m6A <- get_mA_on_ranges(bs_obj, gene_gr = TES)

# Filtrate the TESs that have a coverage of less than 10x
TES_m6A_filt <- TES_m6A %>% filter(mean_coverage >= 10)

# Calculate the mean methylation of the TES
Methylation_TES <- mean(TES_m6A_filt[,6])


#########Unix command line code to produce a repetitive element bed file #########
#
## Build RepeatModeler database using the genome FASTA file
#BuildDatabase -name RepeatModeler_DB.RMdb Genome.fasta
#
## Run RepeatModeler to identify repetitive elements in the genome
#RepeatModeler -pa 6 -LTRStruct -database RepeatModeler_DB.RMdb -genomeSampleSizeMax 100000000
#
## Run RepeatMasker to mask repetitive elements in the genome
#RepeatMasker -nolow -norna -pa 6 -gff -a -xsmall -lib RepeatModeler_DB.RMdb-families.fa Genome.fasta
#
## Extracting relevant columns from RepeatMasker GFF output, formatting, and saving as BED file
#cut Species.RepeatMasker.gff -f 1,4,5,7,9 | cut -d ";" -f1 | awk '{print $1,$2,$3,$5,$5,$4}' | tr " " "\t" > Species.RepeatMasker.bed
#
################################################################################


# Load the .bed file with the repetitive elements
Repeats <- read_bed6_to_GRobject("Cfra.RepeatMasker.bed")

# Calculate the 6mA (regional and per ApT) in the repetitive elements
Repeats_m6A <- get_mA_on_ranges(bs_obj, gene_gr = Repeats)

# Filtrate the repeatitive elements that have a coverage of less than 10x
Repeats_m6A_filt <- Repeats_m6A %>% filter(mean_coverage >= 10)

# Calculate the mean methylation of the repetitive elements
Methylation_Repeats <- mean(Repeats_m6A_filt[,6])


# Calculate the mean methylation value of each  decile group
Methylation_deciles <- aggregate(mApT ~ decile_tpm, data = genes_m6A_filt, FUN = mean)

# Reorder dataframe in descending order
Methylation_deciles <- Methylation_deciles[rev(seq_len(nrow(Methylation_deciles))), ]


# Compute the mean methylation levels of the different regions of all species
Cfra_Multinucleated_Methylation_levels <- matrix(c(name, Methylation_Promoters, Methylation_TSS, Methylation_TES, Methylation_deciles$mApT, Methylation_Repeats))

# Name the columns of the dataframe
Columns_Mean_Methylation <- c("Specie", "Promoters", "TSS+500bs", "TES+250bs", "10th decile", "9th decile", "8th decile", "7th decile", "6th decile", "5th decile", "4th decile", "3rd decile", "2nd decile", "1st decile", "No expression", "Repetitive Elements")

# Create a dataframe with all species
Mean_Methylation_DF <- data.frame(Columns_Mean_Methylation, Cfra_Multinucleated_Methylation_levels, Apar_Methylation_levels, Awhi_Methylation_levels, Cperk_Methylation_levels, Clim_Methylation_levels, Spun_Methylation_levels, Acas_Methylation_levels,
                                Crei_Methylation_levels, Cvar_Methylation_levels, Otau_Methylation_levels, Mpus_Methylation_levels, Alim_Methylation_levels, Ngru_Methylation_levels, Nfow_Methylation_levels, Tvag_Methylation_levels )

# Transpose dataframe
MeanMethylationDF <- transpose(MeanMethylationDF)

# Extract mean methylation levels values
write.csv(MeanMethylationDF, file = 'MeanMeathylationValues.csv', col.names = FALSE,
          row.names = FALSE, sep = "\t")

################################################################################
### AT observed vs expected
################################################################################


##########Calculate the observed versus expected dinucleotides ratios ##########
#
##Run compseq to calculate the number of observed nucleotides in the genome
#compseq Genome.fasta -word 1 Species.1dntps.com
#
##Run compseq to calculate the number of observed dinucleotides in the genome
#compseq Genome.fasta -word 2 Species.2dntps.com
#
##Use the formula (A*T)/Length or (((A+T)/2)^2)/Length to calculate the number of expected dinucleotides in the genome
#
##Divide the number of observed dinucleotides between the number of expected dinucleotides.
#
################################################################################



# Vector with the number of species
Species <- c("Ngru", "Nfow", "Tvag", "Crei", "Cvar", "Ptri", "Otau", "Mpus", "Cpar", "Alim", "Acas", "Ddis", "Spun", "Clim", "Cper", "Awhi", "Cfra", "Aapp")

# Vector with the global percentage of 6mAT of all species
mAT <- c(0.07, 0.11, 3.83, 0.08, 0.17, 0.02, 0.76, 2.57, 0.01, 0.18, 3.47, 0.03, 2.62, 1.24, 4.00, 0.40, 3.01, 1.62)#cat Species.R9_6mA.methyl.context.bed | awk '$17 ~ "[ATCG]AT[ATCG]"' | awk '{cov+=$12; mAT+=$13;} END{print "mAT", cov+mAT, mAT, mAT/(cov+mAT)*100;}' > d

# Vector with the global percentage of 5mCG of all species
mCG <- c(0.43, 0.35, 0.19, 0.30, 30.23, 0.65, 10.51, 5.08, 10.13, 0.19, 0.54, 0.33, 0.44, 0.18, 0.17, 0.74, 0.38, 39.71)#cat Species.R9_5mC.methyl.context.bed | awk '$17 ~ "[ATCG]CG[ATCG]"' | awk '{cov+=$12; mCG+=$13;} END{print "mCG", cov+mCG, mCG, mCG/(cov+mCG)*100;}' > b

# Calculate fraction of the percentage of 5mCG
mCGdivided100 <- mCG/100

# Vector with the values of the observed vs expected ratio of AT
ObservedApTvsExpectedApT <- c(0.959, 0.962, 0.926, 1.000, 1.061, 0.948, 1.082, 1.000, 1.164, 0.892, 0.892, 0.892, 1.012, 0.982, 1.000, 0.897, 1.043, 0.920)

# Vector with the values of the observed vs expected ratio of CG
ObservedCpGvsExpectedCpG <- c(0.423, 0.627, 0.757, 0.895, 0.804,  1.097, 1.545, 1.579, 1.109, 0.861, 0.992, 0.526, 0.899, 0.852, 0.844, 0.846, 0.934, 0.551)


# Dataframe with the percentage od 6mAT and the AT obseverd vs expected ratio
ApTvs6mAT <- data.frame(Species, mAT, ObservedApTvsExpectedApT)

# Dataframe with the percentage od 5mCG and the CG obseverd vs expected ratio
CpGvs5mCG <- data.frame(Species, mCGdivided100, ObservedCpGvsExpectedCpG)

# Prepare dataframe for plot generation
df_longApT <- tidyr::gather(ApTvs6mAT, key = "variable", value = "value", -Species)
df_longCG <- tidyr::gather(CpGvs5mCG, key = "variable", value = "value", -Species)


# Plot bar chart with the 6mAT fraction values and AT observed vs expected ratio
AT_ObservedVsExpected_Methylation_bar_chart <- ggplot(df_longApT, aes(x = reorder(Species, match(Species, Species_order)), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(x = "Species", y = "Values", fill = "Variable") +
  theme_light()
# Plot bar chart with the 5mCG fraction values and CG observed vs expected ratio
CG_ObservedVsExpected_Methylation_bar_chart <- ggplot(df_longCG, aes(x = reorder(Species, match(Species, Species_order)), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(x = "Species", y = "Values", fill = "Variable") +
  theme_light()

Combo_plots_ObservedExpected <- plot_grid(AT_ObservedVsExpected_Methylation_bar_chart, CG_ObservedVsExpected_Methylation_bar_chart, ncol = 1)



################################################################################
### Cell stage specific analysis for A. castellanii and C. fragrantissima
################################################################################

# ------------------------------------------------------------------------------
# Number of methylated genes in different stages
# ------------------------------------------------------------------------------

# Define different cell stages 
Stages <- c("Cfra Multinucleated", "Cfra Amoeba", "Acas Amoeba", "Acas Cyst")


# Vector with the number of genes of the different stages (obtained before with this script)
Genes_stages <- c(8661, 8661, 15529, 15529) 

# Vector with the number of genes with at least 10x coverage of the different stages (obtained before with this script)
Selected_genes_stages <- c(8592, 8355, 15297, 14794) 

# Vector with the number of methylated genes of the different stages with at least 10x coverage (obtained before with this script)
Methylated_genes_stages <- c(7697, 7505, 9717, 9573)

# Calculate the number of unmethylated genes
Unmethylated_genes_stages <- c(Selected_genes_stages- Methylated_genes_stages )

# Calculate the percentage of unmethylated genes
Percentage_of_unmethylated_genes_stages <- c(100*(Unmethylated_genes_stages/Selected_genes_stages))

# Calculate the percentage of methylated genes
Percentage_of_methylated_genes_stages <- c(100*(Methylated_genes_stages/Selected_genes_stages))

# Create a dataframe with the amounts of selected genes, methylated genes and unmethylated genes
Methylated_vs_unmethylated_genes_stages <- data.frame(Stages, Genes_stages, Selected_genes_stages, Methylated_genes_stages, Unmethylated_genes_stages, Percentage_of_unmethylated_genes_stages, Percentage_of_methylated_genes_stages)

#Prepare the dataframes for the plot generation
Percentage_data_long_stages <- gather(Methylated_vs_unmethylated_genes_stages, key = "GeneStatus", value = "Percentage",
                                      Percentage_of_methylated_genes_stages, Percentage_of_unmethylated_genes_stages)
Number_of_genes_data_long_stages <- gather(Methylated_vs_unmethylated_genes_stages, key = "GeneStatus", value = "NumberOfGenes",
                                           Unmethylated_genes_stages, Methylated_genes_stages)


# Plot the percentage of methylation genes
Gene_Methylation_Plot_stages_gg <- ggplot(Percentage_data_long_stages, aes(x = Stages, y = Percentage, fill = reorder(GeneStatus, Percentage))) +
  geom_bar(stat = "identity", position = 'stack') +
  geom_text(aes(label = scales::percent(Percentage, scale = 1)), position = position_stack(vjust = 0.5), color = "black", size = 3) +
  labs(title = "Percentage of Methylated vs Unmethylated Genes",
       x = "Cell Stages",
       y = "Genes Percentage") +
  scale_fill_manual(values = c("#fc8d62", "#66c2a5"), name = "Gene Status",
                    breaks = c('Percentage_of_unmethylated_genes_stages', 'Percentage_of_methylated_genes_stages'),
                    labels = c("Unmethylated Genes", "Methylated Genes")) +
  theme_minimal()

# Plot the number of methylation genes
Number_Gene_Methylation_Plot_stages_gg <- ggplot(Number_of_genes_data_long_stages, aes(x = Stages, y = NumberOfGenes, fill = reorder(GeneStatus, NumberOfGenes))) +
  geom_bar(stat = "identity", position = 'stack') + 
  geom_text(aes(label = NumberOfGenes), position = position_stack(vjust = 0.5), color = "black", size = 3) +
  labs(title = "Number of Methylated vs Unmethylated Genes",
       x = "Cell Stages",
       y = "Number of Genes") +
  scale_fill_manual(values = c("#fc8d62", "#66c2a5"), name = "Gene Status",
                    breaks = c('Unmethylated_genes_stages', 'Methylated_genes_stages'),
                    labels = c("Unmethylated Genes", "Methylated Genes")) +
  theme_minimal()



# ------------------------------------------------------------------------------
# Differentially ApT mehtylated genes in different cell stages
# ------------------------------------------------------------------------------


# Load 6mA methylated genes from different cell stages (obtained before with this script)
Cfrag_Amoeba_genes_6mA <- fread("Cfrag_Amoeba_genes_filt_6mA.tsv")
Cfrag_Multinucleated_genes_6mA <- fread("Cfrag_Multinucleated_genes_filt_6mA.tsv")

# Merge both dataframes
Cfrag_MultVsAmoeba <- merge(Cfrag_Multinucleated_genes_6mA, Cfrag_Amoeba_genes_6mA, by = "gene")

# Calculate difference between methylated ApTs
Cfrag_MultVsAmoeba$delta_6mA_Multi_vs_Amoeba <- Cfrag_MultVsAmoeba$methylated_ApT.x - Cfrag_MultVsAmoeba$methylated_ApT.y 
# Calculate de difference between the methylation fraction 
Cfrag_MultVsAmoeba$delta_6mAfraction_Multi_vs_Amoeba <- Cfrag_MultVsAmoeba$mApT.x - Cfrag_MultVsAmoeba$mApT.y

# Create matrix with the genes and tpm of the different cell stages samples
Cfrag_tpm_matrix <- cbind(fread("abundanceM1.tsv") %>% .[,c(1,5)] %>% dplyr::rename("Multinucleated1" = tpm, gene_id = "target_id") ,
                          fread("abundanceM2.tsv") %>% .[,5] %>% dplyr::rename("Multinucleated2" = tpm),
                          fread("abundanceM3.tsv") %>% .[,5] %>% dplyr::rename("Multinucleated3" = tpm),
                          fread("abundanceA1.tsv") %>% .[,5] %>% dplyr::rename("Amoeba1" = tpm),
                          fread("abundanceA2.tsv") %>% .[,5] %>% dplyr::rename("Amoeba2" = tpm),
                          fread("abundanceA3.tsv") %>% .[,5] %>% dplyr::rename("Amoeba3" = tpm))
# Clean dataframe and remove isoforms
Cfrag_tpm_matrix$gene_id <- substring(Cfrag_tpm_matrix$gene_id, 1, nchar(Cfrag_tpm_matrix$gene_id) - 2)
Cfrag_tpm_matrix <- Cfrag_tpm_matrix[nchar(Cfrag_tpm_matrix$gene_id) <= 10, ]

# check distribution of difference of methylated ApTstpm
Cfrag_MultVsAmoeba$delta_6mA_Multi_vs_Amoeba %>% hist( breaks = 20)

# Plot the histogram of the difference between mApTs in genes
hist(Cfrag_MultVsAmoeba$delta_6mA_Multi_vs_Amoeba, breaks = 20,
     main = "Creolimax    Multinucleated    vs    Amoeba
     Difference Between Genes mApTs",
     xlab = "Difference between mApTs",
     ylab = "Number of Genes")


# Filter the genes that have a difference of 4
Amoeba_hyper_ApTsite <- Cfrag_MultVsAmoeba %>% filter(delta_6mA_Multi_vs_Amoeba < -4) %>% .$gene
Amoeba_hypo_ApTsite <- Cfrag_MultVsAmoeba %>% filter(delta_6mA_Multi_vs_Amoeba > 4) %>% .$gene

# Check number of hyper and hypo methylated genes
length(Amoeba_hyper_ApTsite) 
length(Amoeba_hypo_ApTsite) 

# check distribution of difference of methylated gene body (fraction)
Cfrag_MultVsAmoeba$delta_6mAfraction_Multi_vs_Amoeba %>% hist( breaks = 30)

# Plot the histogram of the difference between the fraction of mApTs in genes
hist(Cfrag_MultVsAmoeba$delta_6mAfraction_Multi_vs_Amoeba, breaks = 30,
     main = "Creolimax      Multinucleated      vs      Amoeba
  Difference Between Genes mApT fraction",
     xlab = "Difference between mApT fraction",
     ylab = "Number of Genes")

# Filter the genes that have a difference of 0.02
Amoeba_hyper_ApTfrac <- Cfrag_MultVsAmoeba %>% filter(delta_6mAfraction_Multi_vs_Amoeba <= -0.02) %>% .$gene
Amoeba_hypo_ApTfrac <- Cfrag_MultVsAmoeba %>% filter(delta_6mAfraction_Multi_vs_Amoeba >= 0.02) %>% .$gene

# Check number of hyper and hypo methylated genes
length(Amoeba_hyper_ApTfrac)
length(Amoeba_hypo_ApTfrac) 

# Check the genes that share status in fraction and in total numbers
table(Amoeba_hyper_ApTsite %in% Amoeba_hyper_ApTfrac)
table(Amoeba_hypo_ApTsite %in% Amoeba_hypo_ApTfrac)

# Add a column to the matrix adding a name to the filtrated genes for clarity  
Cfrag_tpm_matrix_by6mAsite <- Cfrag_tpm_matrix %>% mutate(category = ifelse(gene_id %in% Amoeba_hyper_ApTsite, "Amoeba",
                                                                            ifelse(gene_id %in% Amoeba_hypo_ApTsite, "Multinucleated","Other"))) %>% 
  filter(category != "Other")

# Reshape matrix
Cfrag_tpm_matrix_by6mAsite <- melt(Cfrag_tpm_matrix_by6mAsite)

# Plot the genes selected by the difference of the total amount of mApTs against the methylome of each cell stage
Tpm_by6mAsitediff_gg <- ggplot(Cfrag_tpm_matrix_by6mAsite, aes(x = category, y = log10(value), fill = variable)) + geom_boxplot()


# Same but with the 6mApT fraction
# Add a column to the matrix adding a name to the filtrated genes for clarity
Cfrag_tpm_matrix_by6mAfraction <- Cfrag_tpm_matrix %>% mutate(category = ifelse(gene_id %in% Amoeba_hyper_ApTfrac, "Amoeba",
                                                                                ifelse(gene_id %in% Amoeba_hypo_ApTfrac, "Multinucleated","Other"))) %>% 
  filter(category != "Other")

# Reshape matrix
Cfrag_tpm_matrix_by6mAfraction <- melt(Cfrag_tpm_matrix_by6mAfraction)

# Plot the genes selected by the difference of the fraction of mApTs against the methylome of each cell stage
Cfrag_tpm_matrix_by6mAfraction_gg <- ggplot(Cfrag_tpm_matrix_by6mAfraction, aes(x = category, y = log10(value), fill = variable)) + geom_boxplot()



# ------------------------------------------------------------------------------
# Correlation in the ApT methylation in different cell stages
# ------------------------------------------------------------------------------


# Load the bs object generated before for your first cell stage
Multinucleated_ApT_bs_obj <- readRDS("Cfrag_Multinucleated.AT_bsseq.rds")

# Get the methylation matrix (M) 
Multinucleated_methylated_sites <- getMeth(Multinucleated_ApT_bs_obj, type = "raw", what = "perBase")

# Get the coverage matrix (Cov)
Multinucleated_coverage_matrix <- getCoverage(Multinucleated_ApT_bs_obj)

# Create indices for the high coverage positions
Multinucleated_high_coverage_indices <- which(rowSums(Multinucleated_coverage_matrix) >=  10, arr.ind = TRUE)

# Get the GRanges object from the BSseq object
granges_object <- granges(Multinucleated_ApT_bs_obj)

# Extract the start positions of the methylated sites
Multinucleated_positions <- start(granges_object)
Multinucleated_chromosomes <- seqnames(granges_object)

# Create a dataframe with the percentage of methylation, coverage, and positions
Multinucleated_df <- data.frame(
  Chromosome = Multinucleated_chromosomes[1:length(Multinucleated_positions)],
  Position = Multinucleated_positions,
  Methylation = Multinucleated_methylated_sites,
  Coverage = Multinucleated_coverage_matrix
)

# Change the names of the dataframe columns
names(Multinucleated_df) <- c("Scaffold", "Position", "Methylation", "Coverage")

# Filtrate the sequences that have less than 10x coverage
Multinucleated_df_high_coverage <- Multinucleated_df[Multinucleated_df$Coverage >=  10, ]

# Load the bs object generated before for your second cell stage
Amoeba_ApT_bs_obj <- readRDS("../Cfrag_Amoeba.AT_bsseq.rds")

# Get the methylation matrix (M)
Amoeba_methylated_sites <- getMeth(Amoeba_ApT_bs_obj, type = "raw", what = "perBase")

# Get the coverage matrix (Cov)
Amoeba_coverage_matrix <- getCoverage(Amoeba_ApT_bs_obj)

# Create indices for the high coverage positions
Amoeba_high_coverage_indices <- which(rowSums(Amoeba_coverage_matrix) >=  10, arr.ind = TRUE)

# Get the GRanges object from the BSseq object
Amoeba_granges_object <- granges(Amoeba_ApT_bs_obj)

# Extract the start positions of the methylated sites
Amoeba_positions <- start(Amoeba_granges_object)
Amoeba_chromosomes <- seqnames(Amoeba_granges_object)


# Create a dataframe with the percentage of methylation, coverage, and positions
Amoeba_df <- data.frame(
  Chromosome = Amoeba_chromosomes[1:length(Amoeba_positions)],
  Position = Amoeba_positions,
  Methylation = Amoeba_methylated_sites,
  Coverage = Amoeba_coverage_matrix
)

# Change the names of the dataframe columns
names(Amoeba_df) <- c("Scaffold", "Position", "Methylation", "Coverage")

# Filtrate the sequences that have less than 10x coverage
Amoeba_df_high_coverage <- Amoeba_df[Amoeba_df$Coverage >=  10, ]


# Merge the data frames
Amoeba_Multinucleated_df <- merge(Amoeba_df_high_coverage, Multinucleated_df_high_coverage, by = c("Scaffold", "Position"), all = TRUE)

# Replace NA values with 0
Amoeba_Multinucleated_df[is.na(Amoeba_Multinucleated_df)] <- 0

# Rename columns for clarity
names(Amoeba_Multinucleated_df)[names(Amoeba_Multinucleated_df) == "Methylation.x"] <- "Amoeba_Methylation"
names(Amoeba_Multinucleated_df)[names(Amoeba_Multinucleated_df) == "Methylation.y"] <- "Multinucleated_Methylation"

# Merge the scaffold and nucleotide position
Amoeba_Multinucleated_df$Position <- paste(Amoeba_Multinucleated_df$Scaffold, ":", Amoeba_Multinucleated_df$Position)

# Filtrate the sequences that have less than 10x coverage in any of the samples
Amoeba_Multinucleated_df_high_coverage <- Amoeba_Multinucleated_df[Amoeba_Multinucleated_df$Coverage.x >=  10 & Amoeba_Multinucleated_df$Coverage.y >=  10, ]

# Smooth scatter the ApT methylation from the different cell stages
smoothScatter(Amoeba_Multinucleated_df_high_coverage$Amoeba_Methylation,Amoeba_Multinucleated_df_high_coverage$Multinucleated_Methylation,
              xlab = "Amoeba ApTs methylation percentage per site", ylab = "Multinucleated ApTs methylation percentage per site")
cor(Amoeba_Multinucleated_df_high_coverage$Amoeba_Methylation,Amoeba_Multinucleated_df_high_coverage$Multinucleated_Methylation, method = "pearson") 




