# Intall the required libraries

install.packages("data.table")
install.packages("stringr")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("tidyr")
install.packages("R.utils")
install.packages("reshape2")
install.packages("ggseqlogo")
install.packages("cowplot")
install.packages("BiocManager")
BiocManager::install("bsseq")

# Load the required libraries

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
Negative_species <- c("Ddis", "Cpar", "Ptri")

# ------------------------------------------------------------------------------
# Load genomic information and filtrate file
# ------------------------------------------------------------------------------


#########Unix command line code to produce a genomic file from nanopore #########
#
#########For data generated with R9 flowcells
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
#bedtools getfasta -fi Genome.fasta -fo Species.R9_6mA.methyl.plus3.fa -s -bed Species.R9_6mA.methyl.plus3.bed
#
## Run seqkit to convert FASTA format to tab-separated format and extracting sequences
#seqkit fx2tab Species.R9_6mA.methyl.plus3.fa | cut -f 2 > a
#
## Merge BED file containing genomic coordinates with the extracted sequences
#paste Species.R9_6mA.methyl.bed a > Species.R9_6mA.methyl.context.bed
#
## Filter lines containing the pattern [ATCG]AT[ATCG] and extract columns 1, 2, 3, and 11 to produce a bedgraph file
#awk '$17 ~ /[ATCG]AT[ATCG]/' Species.R9_6mA.methyl.context.bed | cut -f 1,2,3,11 > Species.R9_6mA.methyl.ApT.bedGraph
#
## Sort the bedGraph file
#LC_COLLATE=C sort -k1,1 -k2,2n Species.R9_6mA.methyl.ApT.bedGraph > Species.R9_6mA.methyl.ApT.sorted.bedGraph
#
## Run bedGrpahToBigWig to convert the sorted bedGraph file into bigwig
#bedGraphToBigWig Species.R9_6mA.methyl.ApT.sorted.bedGraph Genome.fasta.fai Species.R9_6mA.methyl.ApT.bigwig
#
## Bigwig files are used later on the genome browser (IGV) and in DeepTools2 plots
#
#
#########For data generated with R10 flowcells
#
## Run dorado to do the basecalling 
#
#dorado basecaller sup,6mA pod5/ --reference Genome.fasta > Allmods.bam
#
## Run samtools to sort the BAM file
#samtools sort -@ 20 Allmods.bam -o Allmods.sorted.bam
#
## Run samtools to index the sorted BAM file
#samtools index -@ 20 Allmods.sorted.bam
#
## Run modkit to convert the sorted BAM file to BED format
#
#modkit pileup --ref Genome.fasta --motif AN 0 --mod-thresholds a:0.995 Allmods.sorted.bam Species.R10.pileup.AN.0995.bed --log-filepath Species.R10.pileup.AN.log
#modkit pileup --ref Genome.fasta --motif AT 0 --mod-thresholds a:0.995 Allmods.sorted.bam Species.R10.pileup.AT.0995.bed --log-filepath Species.R10.pileup.AT.log
#
## Extract columns 1, 2, 3, and 11 to produce a bedgraph file
#awk -v OFS='\t' '{print $1, $2, $3, $11}'  Species.R10.pileup.AN.0995.bed > Cpecies.R10_6mA.methyl.ApT.0995.bedGraph
#
## Run bedGrpahToBigWig to convert the sorted bedGraph file into bigwig
#bedGraphToBigWig Cpecies.R10_6mA.methyl.ApT.0995.bedGraph Genome.fasta.fai Species.R10_6mA.methyl.ApT.bigwig
#
## Bigwig files are used later on the genome browser (IGV) and in DeepTools2 plots
#
################################################################################



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


# ------------------------------------------------------------------------------
# Plot the correlation values of all species
# ------------------------------------------------------------------------------


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


# Generate the bar chart of correlation values of each species
AT_Pearson_Correlation_Bar_Chart<- ggplot(AT_Pearson_Correlation_DF_long, aes(x = reorder(Species, match(Species, Species_order)), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(title = "ApT symmetric methylation correlation", x = 'Species', y = "Pearson correlation value", fill = "Variable") +
  theme_light() +
  guides(fill = "none") +
  scale_y_continuous(limits = c(min(df_longA$value), max(df_longA$value)))



# ------------------------------------------------------------------------------
# Load gene bed file and get the 6mA from it
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

# Calculate the 6mA (regional and per ApT) for every gene
genes_m6A <- get_mA_on_ranges(bs_obj = bs_obj_collapsed, gene_gr = genes)


# Check the coverage of each gene 
genes_m6A$mean_coverage %>% hist( breaks = 50, main = "Mean coverage per gene")

# Discard genes with mean coverage < 10x 
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

# Create a dataframe for each species with the filtrated 6mA dataset
Cfra_df <- genes_m6A_filt 

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

# Plot the percentage of methylation genes
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


# Plot the numnber of methylation genes
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
# Calculate number genes per methylation category 
# ------------------------------------------------------------------------------

# Combine the 6mA filtered data frames of each species into a single one
combined_genes_m6A_filt <- bind_rows(Otau_df, Acas_df, Mpus_df, Alim_df, Aapp_df, Awhi_df, Cvar_df, Cfra_df, Clim_df, Cper_df, Crei_df, Nfow_df, Ngru_df, Spun_df, Tvag_df, Ddis_df, Ptri_df, Cpar_df)

# Separate genes according to different methylation levels
genes_by_m6ApT_levels_groups <- combined_genes_m6A_filt %>%
  mutate(mApT_bin = case_when(
    mApT <= 0.01 ~ "0 - 1%",
    mApT <= 0.05 ~ "1 - 5%",
    mApT <= 0.10 ~ "5 - 10%",
    TRUE        ~ "> 10%"
  )) %>%
  group_by(mApT_bin, Species) %>%
  summarise(gene_count = n(), .groups = "drop")

# Establish the order of the groups for plotting
genes_by_m6ApT_levels_groups <- genes_by_m6ApT_levels_groups %>%
  mutate(mApT_bin = factor(mApT_bin, levels = c("> 10%", "5 - 10%", "1 - 5%", "0 - 1%")))


# Plot the fraction of genes per 6mApT level
Gene_number_per_mAPT_levels_fraction <- ggplot(genes_by_m6ApT_levels_groups, aes(x = reorder(Species, match(Species, All_species_order)), y = gene_count, fill = mApT_bin)) +
  geom_bar(position="fill", stat="identity") +
  labs(x = "Species", y = "Fraction of Genes", title = "Fraction of genes per global mApT% levels", fill = "Gene mApT% levels") +
  #scale_fill_grey(start = 0.7, end = 0.1) + # Use grey scale
  scale_fill_manual(
    values = c("0 - 1%" = "#105965FF", "1 - 5%" = "#217A79FF", "5 - 10%" = "#6CC08BFF", "> 10%" = "#D3F2A3FF")  # Specify colors for each category
  ) +
  theme_minimal()


# ------------------------------------------------------------------------------
# Calculate number and percentage of methylated As and ATs
# ------------------------------------------------------------------------------


# Create a vector with the number of adenines with a coverage >=10x of all species
Adenine10cov <- c(23365911, 90853275, 10895968, 18220779, 109724, 11500215, 17433674, 9385095, 37942246, 68793, 5005691, 170572, 2380137, 13856641, 32187353, 23423770, 18367077, 65745344)#awk '$10 >= 10' Ngru.R9_6mA.methyl.context.bed > Ngru.AT.Cov.bed

# Create a vector with the number of adenines with a coverage >=10x of all species and that are methylated at lest 10% of the times of all species
MethyAdenine10cov10perc <- c(944769, 2827469, 186850, 1059771, 2218, 380685, 454255, 192603, 22183, 133, 56598, 4642, 1670, 221, 72974, 21894, 13737, 2632920)#awk '$10 >= 10 && $11 >= 10' Ngru.R9_6mA.methyl.context.bed > Ngru.AT.CovMeth.bed

# Create a vector with the number of ATs with a coverage >=10x of all species
AT10cov <- c(7185587, 20046251, 3665666, 4662423, 28652, 2827680, 3217753, 3242009, 5332264, 13571, 888790, 21886, 485880, 3051863, 7884208, 7518561, 4244237, 11732126)#awk '$NF ~ /[ATCG]AT[ATCG]/ && $10 >= 10' Ngru.R9_6mA.methyl.context.bed > Ngru.AT.Cov.bed

# Create a vector with the number of ATs with a coverage >=10x of all species and that are methylated at lest 10% of the times of all species
MethAT10cov10perc <- c(554769, 1366197, 64270, 632364, 856, 209721, 319880, 11333, 15954, 66, 24988, 1790, 230, 29, 45515, 14533, 9669, 1233568)#awk '$NF ~ /[ATCG]AT[ATCG]/ && $10 >= 10 && $11 >= 10' Ngru.R9_6mA.methyl.context.bed > Ngru.AT.CovMeth.bed

# Calculate the percentage of methylated As
PercentageMethylA <- c(100*(MethyAdenine10cov10perc/Adenine10cov))

# Calculate the percentage of methylated ATs
PercentageMethylAT <- c(100*(MethAT10cov10perc/AT10cov))

# Create a dataframe with this values
PercentageMethylatedAdenines <- data.frame(Species_order, Adenine10cov, MethyAdenine10cov10perc, AT10cov, MethAT10cov10perc, PercentageMethylA, PercentageMethylAT)


#Create bar plots with the percentage of methylated A and ATs
Percentage_Methy_A_bar_chart <- ggplot(PercentageMethylatedAdenines, aes(x = Species_order, y = PercentageMethylA, fill = "#66c2a5")) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "Percentage of Methylated Adenines", x= 'Species', y = "%") +
  geom_text(aes(label = Adenine10cov), position = position_stack(vjust = 0.5), color = "black", size = 2) +
  scale_fill_manual(values = "#66c2a5")  +
  theme_light() +
  guides(fill = "none")
Percentage_Methy_AT_bar_chart <- ggplot(PercentageMethylatedAdenines, aes(x = Species_order, y = PercentageMethylAT, fill = "#66c2a5")) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "Percentage of Methylated ApT", x= 'Species', y = "%") +
  geom_text(aes(label = AT10cov), position = position_stack(vjust = 0.5), color = "black", size = 2) +
  scale_fill_manual(values = "#66c2a5")  +
  theme_light() +
  guides(fill = "none")
combo_plots <- plot_grid(Percentage_Methy_A_bar_chart, Percentage_Methy_AT_bar_chart, ncol = 1 )


#Prepare the dataframe for the plot generation
Number_of_methylated_sites_data_long <- gather(PercentageMethylatedAdenines, key = "SiteStatus", value = "NumberOfSites",
                                               Adenine10cov, MethyAdenine10cov10perc)
Number_of_methylated_ATsites_data_long <- gather(PercentageMethylatedAdenines, key = "SiteStatus", value = "NumberOfSites",
                                                 AT10cov, MethAT10cov10perc)

#Create bar plots with the number of methylated A and ATs
Number_Methy_A_bar_chart <- ggplot(Number_of_methylated_sites_data_long, aes(x = Species_order, y = NumberOfSites, fill = reorder(SiteStatus, -NumberOfSites))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "Number of Adenines and Methylated Adenines", x = 'Species', y = "log(Number of Adenenines)") +
  scale_fill_manual(values = c("#fc8d62", "#66c2a5"), name = "Site Status",
                    breaks = c('Adenine10cov', 'MethyAdenine10cov10perc'),
                    labels = c("Number of Adenines", "Number of Methylated Adenines")) +
  theme_light() +
  scale_y_log10()
Number_Methy_AT_bar_chart <- ggplot(Number_of_methylated_ATsites_data_long, aes(x = Species_order, y = NumberOfSites, fill = reorder(SiteStatus, -NumberOfSites))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "Number of ApTs and Methylated ApTs", x = 'Species', y = "log(Number of ApTs)") +
  scale_fill_manual(values = c("#fc8d62", "#66c2a5"), name = "Site Status",
                    breaks = c('AT10cov', 'MethAT10cov10perc'),
                    labels = c("Number of ApTs", "Number of Methylated ApTs")) +
  theme_light()  +
  scale_y_log10()
combo_plots <- plot_grid(bar_chartANumber, bar_chartATNumber, ncol = 1 )


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
#cat Species.R9_6mA.methyl.context.bed | awk '$17 ~ "TAT[CGT]"' | awk '{cov+=$12; TmAT+=$13;} END{print "TmAT", cov+TmAT, TmAT, TmAT/(cov+TmAT)*100;}' > k
#
# Create a header for the output file
#
#echo "context coverage m m_perc" > aa
#
# Concatenate all result files into one file
#
#cat aa a b c d e f g h i j k > Species_R9_methyl_stats.txt
#
################################################################################

# Load 6mA percentages levels calculated before
df <- fread("Species_R9_methyl_stats.txt")

# Organize data frame according to species
df$Species <- factor(df$Species, levels = df$Species)

# Select only the dinucleotide context of each specie
Dinucleotide_6mA_percentage <- df %>%
  pivot_longer(cols = c(mAT, mAA, mAC, mAG), names_to = "Context", values_to = "Percentage")

# Select the rest of the context  of each specie
Extended_context_6mA_percentage <- df %>%
  pivot_longer(cols = c(SmATS, WmATW, VmATB, TmAT), names_to = "Context", values_to = "Percentage" )

# Establish a context order for the plots
Dinucleotide_context_order <- c("mAT", "mAA", "mAG", "mAC")
Tetranucleotide_context_order <- c("SmATS", "VmATB", "WmATW", "TmAT")

# Plot bar chart of 6mA percentages
mALevels <- ggplot(data = df, aes(x = reorder(Species, match(Species, All_species_order)), y = mA)) +
  geom_bar(stat = "identity", fill="darkgrey") +
  labs(y = "percentage 6mA (%)") +
  theme_light()

# Plot stacked bar chart of 6mA dinucleotide fraction levels
Fraction_methylated_dinucleotide_stacked_bar_chart <- ggplot(Dinucleotide_6mA_percentage, aes(x = reorder(Species, match(Species, All_species_order)), y = Percentage, fill = Context)) +
  geom_bar(position="fill", stat="identity") +
  labs(title = "Fraction of Methylated Dinucleotides on Total 6mA ",
       y = "Fraction of Total 6mA ", x = "Species") +
  scale_fill_manual(values = c ("#217A79FF", "#D3F2A3FF", "#6CC08BFF", "#074050FF")) +
  theme_minimal()

# Create bar chart of the 6mA dinucleotide percentages 
Dinucleotide_6mA_percentage_bar_chart <- ggplot(Dinucleotide_6mA_percentage, aes(x = reorder(Species, match(Species, All_species_order)), y = Percentage, fill = reorder(Context, match(Context, Dinucleotide_context_order)))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) + 
  labs(title = "Dinucleotide 6mA %",
       y = "6mA % per Dinucleotide", x = "Species", fill = "Context") +
  scale_fill_manual(values = c("#074050FF", "#217A79FF", "#6CC08BFF", "#D3F2A3FF")) +
  theme_minimal()

# Create bar chart zooming in the 6mA dinucleotide percentages of the negative species
Dinucleotide_zoom_bar_chart <- ggplot(filter(Dinucleotide_6mA_percentage, Species %in% Negative_species), 
                                aes(x = reorder(Species, match(Species, All_species_order)), 
                                    y = Percentage, fill = reorder(Context, match(Context, Dinucleotide_context_order)))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  coord_cartesian(ylim = c(0, max(filter(Dinucleotide_6mA_percentage, Species %in% Negative_species)$Percentage))) + 
  labs(title = "Negative Species Zoomed-in View", x = NULL, y = NULL) +
  scale_fill_manual(values = c("#074050FF", "#217A79FF", "#6CC08BFF", "#D3F2A3FF")) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", color = NA),   # White background
        plot.background = element_rect(fill = "white", color = NA))    # White plot background

# Combine both plots
Dinucleotide_context_6mA_percentage_and_zoom_bar_chart <- ggdraw() +
  draw_plot(Dinucleotide_6mA_percentage_bar_chart) +
  draw_plot(Dinucleotide_zoom_bar_chart, x = 0.5, y = 0.7, width = 0.3, height = 0.3) 

# Repeat to plot the extended context
Extended_context_6mA_percentage_bar_chart <- ggplot(Extended_context_6mA_percentage, aes(x = reorder(Species, match(Species, All_species_order)), y = Percentage, fill = reorder(Context, match(Context, Tetranucleotide_context_order)))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) + 
  labs(title = "Extended Context 6mA %",
       y = "6mA % per Sequence", x = "Species", fill = "Context") +
  scale_fill_manual(values = c("#D82632FF", "#F76D5EFF", "#FFAD72FF", "#FFE099FF")) +
  theme_minimal()

Extended_zoom_bar_chart <- ggplot(filter(Extended_context_6mA_percentage, Species %in% Negative_species), 
                   aes(x = reorder(Species, match(Species, All_species_order)), 
                       y = Percentage, fill = reorder(Context, match(Context, Tetranucleotide_context_order)))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  coord_cartesian(ylim = c(0, max(filter(Extended_context_6mA_percentage, Species %in% Negative_species)$Percentage))) + 
  labs(title = "Negative Species Zoomed-in View", x = NULL, y = NULL) +
  scale_fill_manual(values = c("#D82632FF", "#F76D5EFF", "#FFAD72FF", "#FFE099FF")) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", color = NA),   
        plot.background = element_rect(fill = "white", color = NA))    

Extended_context_6mA_percentage_and_zoom_bar_chart <- ggdraw() +
  draw_plot(Extended_context_6mA_percentage_bar_chart) +
  draw_plot(Extended_zoom_bar_chart, x = 0.6, y = 0.6, width = 0.3, height = 0.3) 



################################################################################
###Comparison between A methylation and transcription
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
##Run stringtie to assemble transcripts 
#stringtie --rf -G Genomic.gtf -A Species.RNAseq_abundance.tsv -e ${name}.hisat2.sorted.bam
#
################################################################################



# Load RNA-seq data
tpm <-fread("Cfra_Multinucleated.RNAseq_abundance.tsv") %>% dplyr::rename(gene = "Gene\ ID") %>% 
  dplyr::select(gene, TPM)

# Loin RNA and 6mA data
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
# These files were later plotted with DeepTools2 computeMatrix function

####Unix command line code to plot different gene expression deciles against 6mAT levels####
#
### Extract the first column of decile files generated before
#awk '{print $1}' Species0-10%.txt > first_decile
#awk '{print $1}' Species10-20%.txt > second_decile
#…Repeat for the rest of deciles…
#awk '{print $1}' Species90-100%.txt > tenth_decile
#awk '{print $1}' SpeciesNoexpr.txt > no_expression
#
## Extract gene annotations corresponding to the genes in each decile
#fgrep -w -f first_decile Species.annotations.gene.bed > first_decile.bed
#fgrep -w -f second_decile Species.annotations.gene.bed > second_decile.bed
#…repeat for the rest of deciles…
#fgrep -w -f tenth_decile Species.annotations.gene.bed > tenth_decile.bed
#fgrep -w -f no_expression Species.annotations.gene.bed > no_expression.bed
#
## Compute matrix using the Transcription Start Site as reference point
#computeMatrix reference-point --referencePoint TSS -out Species.TSS.expression.dtmatrix -S Species.R9_6mA.methyl.ApT.bigwig -R tenth_decile.bed seventh_decile.bed sixth_decile.bed third_decile.bed second_decile.bed no_expression.bed -a 2000 -b 2000 -bs 10 -p 20
#
## Plot heatmap using the computed matrix
#plotHeatmap -m Species.TSS.expression.dtmatrix -out Species.TSS.expression.heatmap.pdf --colorMap Reds --regionsLabel 90-100% 60-70% 50-60% 20-30% 10-20% "No expression" --samplesLabel 6mAT --missingDataColor="silver" --interpolationMethod nearest
#
################################################################################

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


# ------------------------------------------------------------------------------
# Number of methylated genes by decile of expression
# ------------------------------------------------------------------------------

# Asign methylation status to genes
genes_m6A_filt <- genes_m6A_filt %>% mutate(status = ifelse(methylated_ApT >= 3, "Methylated","Unmethylated"))

# Count number of methylated and unmethylated genes per decile of expression
methylation_counts <- genes_m6A_filt %>%
  group_by(decile_tpm, status) %>%
  summarise(count = n(),
            mean_mApT = ifelse(status == "Methylated", mean(mApT, na.rm = TRUE), NA)) %>%
  ungroup()
methylation_counts$status[is.na(methylation_counts$status)] <- "Unmethylated"

# Plot the number of methylated genes per decile of expression and assign methylation levels
Number_and_Levels_of_Methylated_Genes <- ggplot(methylation_counts,
                                                aes(x = factor(decile_tpm), y = count, group = status)) +
  geom_bar(aes(fill = mean_mApT),
           stat = "identity",
           position = position_dodge(width = 0.6),
           alpha = 1) +
  scale_fill_gradient(low = "#fdf9f7", high = "#8a1811") +
  labs(title = "Gene Methylation Status by Expression Decile in C.fra",
       x = "Expression Decile",
       y = "Number of Genes",
       fill = "Mean mApT Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



################################################################################
###Check promoter methylation
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
#these files were later plotted with DeepTools2 computeMatrix function

##Unix command line code to plot Methylated and Unmethylated genes bed files against 6mAT levels##
#
## Compute matrix using the Transcription Start Site as reference point from Deeptools2
#computeMatrix reference-point --referencePoint TSS -out SpeciesTSS.MethylatedGenes.dtmatrix -S Species.R9_6mA.methyl.ApT.bigwig -R Species.Methyl.bed Species.Unmethyl.bed -a 2000 -b 2000 -bs 10 -p 20
#
## Plot heatmap using the computed matrix
#plotHeatmap -m SpeciesTSS.MethylatedGenes.dtmatrix -out SpeciesTSS.MethylatedGenes.heatmap.pdf --colorMap Reds --regionsLabel "m6A Genes" "other Genes" --samplesLabel 6mAT --missingDataColor="silver" --interpolationMethod nearest
#
################################################################################ 


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
#these files were later plotted with DeepTools2 computeMatrix function


##Unix command line code to plot gene orientation bed files against 6mAT levels#
#
## Compute matrix using the Transcription Start Site as reference point
#computeMatrix reference-point --referencePoint TSS -out SpeciesTSS.GeneOrientation.dtmatrix -S Species.R9_6mA.methyl.ApT.bigwig -R Species.head_to_head.methyl.bed Species.head_to_tail.methyl.bed -a 2000 -b 2000 -bs 10 -p 20
#
## Plot heatmap using the computed matrix
#plotHeatmap -m SpeciesTSS.GeneOrientation.dtmatrix -out SpeciesTSS.GeneOrientation.heatmap.pdf --colorMap Reds --regionsLabel “Head to head methylated genes" “Head to tail methylated genes" --samplesLabel 6mAT --missingDataColor="silver" --interpolationMethod nearest
#
################################################################################ 



################################################################################
###CG correlation plot
################################################################################

# Load the .bed file with the genomic information generated before
mC_df <- fread("../Cfrag.R9_5mC.methyl.context.bed.gz")

# Get the % of CG
CG <- mC_df %>% filter(grepl(V17, pattern = "[ATCG]CG[ATCG]"))

# Extract CG information from genomic data
CG_dat <- CG %>% dplyr::rename(chr = V1, start = V2, end = V3, strand = V6, mC_perc = V11,
                                C_reads = V13) %>% mutate(CT_reads = V12 + C_reads)

# Calculate correlation between + and - strands 
CG_dinuc <- CG_dat %>% mutate(locus = if_else(strand == "+", paste0(chr,":",start+1),paste0(chr,":",start) )) %>%
  dplyr::select(locus, strand, mC_perc, CT_reads)

# Filter positions with 10x coverage in both strands
CG_dinuc_watson <- CpG_dinuc %>% filter(CT_reads >= 10, strand == "+")
CG_dinuc_crick <- CpG_dinuc %>% filter(CT_reads >= 10, strand == "-")

# Take only those present in both data frames (10x coverage in both strands)
CG_dinuc <- inner_join(CG_dinuc_watson, CG_dinuc_crick, by = "locus")

# Calculate Pearson and Spearman correlation

CG_pearson <- cor(CG_dinuc$mC_perc.x,CG_dinuc$mC_perc.y, method = "pearson")

CG_spearman <- cor(CG_dinuc$mC_perc.x,CG_dinuc$mC_perc.y, method = "spearman")

# Compute the correlation values for all species
Cfra_CGCorrelation <- c(name, pearson, spearman )#Repeat with the rest of species


# Prepare a matrix with all species values
CG_Pearson_Correlation <- matrix(
  c(
    Cfra_CGCorrelation[1], Cfra_CGCorrelation[2],
    Aapp_CGCorrelation[1], Aapp_CGCorrelation[2],
    # ... Repeat for the rest of species
    Tvag_CGCorrelation[1], Tvag_CGCorrelation[2]
  ),
  ncol = 2,
  byrow = TRUE
)

# Prepare the dataframe for plot generation
CG_Pearson_Correlation_DF <- as.data.frame(CG_Pearson_Correlation)
colnames(CG_Pearson_Correlation_DF) <- c("Species", "Pearson")
CG_Pearson_Correlation_DF_long <- tidyr::gather(CG_Pearson_Correlation_DF, key = "variable", value = "value", -Species)
CG_Pearson_Correlation_DF_long$value <- as.numeric(CG_Pearson_Correlation_DF_long$value)

# Generate the bar chart of the correlation values of each species
CG_Pearson_Correlation_Bar_Chart<- ggplot(CG_Pearson_Correlation_DF_long, aes(x = reorder(Species, match(Species, Species_order)), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(title = "CpG symmetric methylation correlation", x = 'Species', y = "Pearson correlation value", fill = "Variable") +
  theme_light() +
  ylim(-0.04,0.70) +
  guides(fill = "none") +
  scale_y_continuous(limits = c(min(df_longA$value), max(df_longA$value)))




################################################################################
###Mean methylation levels across context
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


# Load the .bed file with the repeatiteve elements
Repeats <- read_bed6_to_GRobject("Cfra.RepeatMasker.bed")

# Calculate the 6mA (regional and per ApT) in the repeatitive elements
Repeats_m6A <- get_mA_on_ranges(bs_obj, gene_gr = Repeats)

# Filtrate the repeatitive elements that have a coverage of less than 10x
Repeats_m6A_filt <- Repeats_m6A %>% filter(mean_coverage >= 10)

# Calculate the mean methylation of the repeatitive elements
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
###AT and CG observed vs expected
################################################################################


##########Calculate the observed versus expected dinucleotides ratios ##########
#
## Run compseq to calculate the number of observed nucleotides in the genome
#compseq Genome.fasta -word 1 Species.1dntps.com
#
## Run compseq to calculate the number of observed dinucleotides in the genome
#compseq Genome.fasta -word 2 Species.2dntps.com
#
## Use the formula (A*T)/Length or (((A+T)/2)^2)/Length to calculate the number of expected dinucleotides in the genome
#
## Divide the number of observed dinucleotides between the number of expected dinucleotides.
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


# Dataframe with the percentage of 6mAT and the AT obseverd vs expected ratio
ApTvs6mAT <- data.frame(Species, mAT, ObservedApTvsExpectedApT)

# Dataframe with the percentage of 5mCG and the CG obseverd vs expected ratio
CpGvs5mCG <- data.frame(Species, mCGdivided100, ObservedCpGvsExpectedCpG)

# Prepare dataframe for plot generation
df_longApT <- tidyr::gather(ApTvs6mAT, key = "variable", value = "value", -Species)
df_longCG <- tidyr::gather(CpGvs5mCG, key = "variable", value = "value", -Species)


# Plot bar chart of the 6mAT fraction values and AT observed vs expected ratio
AT_ObservedVsExpected_Methylation_bar_chart <- ggplot(df_longApT, aes(x = reorder(Species, match(Species, Species_order)), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(x = "Species", y = "Values", fill = "Variable") +
  theme_light()
# Plot bar chart of the 5mCG fraction values and CG observed vs expected ratio
CG_ObservedVsExpected_Methylation_bar_chart <- ggplot(df_longCG, aes(x = reorder(Species, match(Species, Species_order)), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(x = "Species", y = "Values", fill = "Variable") +
  theme_light()

Combo_plots_ObservedExpected <- plot_grid(AT_ObservedVsExpected_Methylation_bar_chart, CG_ObservedVsExpected_Methylation_bar_chart, ncol = 1)


################################################################################
###Genome composition analysis
################################################################################

##Unix command line code to plot AT, TAT and ATA  in Methylated and Unmethylated genes bed files##
#
## Run megalodon_extras to create BED file containing motif positions for motif "AT" (Change according to the motif of interest --motif TAT --motif ATA )
#megalodon_extras modified_bases create_motif_bed --motif AT 0 --out-filename Cfra.AT-motif.consensus.bed Genome.fasta
#
## Convert the BED file to bedGraph format, adjusting coordinates for positive strand
#cat Cfra.AT-motif.consensus.bed | awk '$6 == "+" {print $1,$2,$3+1,"1"}' > Cfra.AT-motif.consensus.bedGraph
#
## Sort the bedgraph file
#LC_COLLATE=C sort -k1,1 -k2,2n Cfra.AT-motif.consensus.bedGraph > Cfra.AT-motif.consensus.sorted.bedGraph
#
## Run bedGrpahToBigWig to convert the sorted bedGraph file into bigwig
#bedGraphToBigWig Cfra.AT-motif.consensus.sorted.bedGraph Genome.fasta.fai Cfra.AT-motif.consensus.bigwig
#
## Use bigwigs files to plot against methylated and unmethylated genes using Deeptools2 as done before
#
#computeMatrix reference-point --referencePoint TSS -out Cfra.AT.TSS.composition.dtmatrix \
#-S \
#Cfra.AT-motif.consensus.bigwig \
#-R Cfra.Methyl.bed Cfra.Unmethyl.bed \
#-a 2000 -b 2000 -bs 10 -p 20 --missingDataAsZero
#
# Plot heatmap using the computed matrix
#
#plotHeatmap -m Cfra.AT.TSS.composition.dtmatrix -out Cfra.AT.TSS.composition.heatmap.pdf \
#--colorMap Reds \
#--samplesLabel AT \
#--missingDataColor="silver" --interpolationMethod nearest
#
################################################################################


################################################################################
###Stage-specific transcriptional changes vs 6mA analysis
################################################################################


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

Cfra_AmoebavsMultinucleated_Pearson <- cor(Amoeba_Multinucleated_df_high_coverage$Amoeba_Methylation,Amoeba_Multinucleated_df_high_coverage$Multinucleated_Methylation, method = "pearson")


# ------------------------------------------------------------------------------
# Differentially expressed genes in different cell stages
# ------------------------------------------------------------------------------



#########Unix command line code to produce an abundance transcriptomic file with counts from raw RNA-seq fastq file#########
#
## Run gffread to create a mRNA.fasta file
#gffread -w Species.mRNA.fasta -g Genome.fasta Species.annotation.gff3
#
## Run kallisto to generate index ile
#kallisto index --index=Species.index Species.mRNA.fasta
#
## Run kallisto to create abundance files with counts numbers
#kallisto quant -i Species.index -o Species.counts.tsv Species_RNA_seq_1.fastq Species_RNA_seq_2.fastq
#
################################################################################


# Read and reformat RNA-seq file with counts to make a matrix
Counts <- fread("CfragMultinucleatedAmoeba.counts.tsv")
Counts <- tpm[grep(Counts$target_id, pattern = "\\.1"),]

# Remove de different isoforms choosing the elements that are longer than 10 characters (This will change depending on the annotation of each species)
Counts <- Counts[nchar(Counts$target_id) <= 10, ] 
Counts$gene_id <- substring(Counts$target_id, 1, nchar(Counts$target_id) - 2)

Counts_mtx <- Counts %>% dplyr::select(-gene_id, -target_id)
rownames(Counts_mtx) <- Counts$gene_id

# Filter genes that are not expressed
allZero <- rowSums(Counts_mtx==0)==ncol(Counts_mtx)
gene_counts <- Counts_mtx[!allZero,]
rownames(gene_counts ) <- rownames(Counts_mtx)[!allZero] %>% as.character()
gene_counts <- round(gene_counts, digits = 0)

# Plot the sample correlation to see if replicates make sense
d.correlation <- as.dist(1 - cor(gene_counts,method=c("spearman")))
fit <- hclust(d.correlation, method="complete")
pdf("Cfrag_RNAseq_clustering_byDEgenes.pdf", height = 5)
plot(fit) # display dendogram
dev.off()


# Generate the experimental conditions data frame
coldata <- data.frame(sample=colnames(gene_counts)) %>%
  mutate(condition = ifelse(grepl("Multinucleated", sample), "Multinucleated",
                            ifelse(grepl("Amoeba", sample), "Amoeba","Broad"))) 

# Create the DEseq2 object
Dds_all <- DESeqDataSetFromMatrix(countData = gene_counts,
                                  colData = coldata,
                                  design = ~ condition)
rownames(Dds_all) <- rownames(gene_counts)
Dds_all <- estimateSizeFactors(Dds_all)

# Filter very lowly expressed genes
Dds_all <- Dds_all[rowSums(Counts(Dds_all)) >= 10,]
Dds_all <- Dds_all[rowSums(Counts(Dds_all)==0) <= ncol(Counts(Dds_all))*0.5,]

# Function to get differentially expressed genes
Get_diff_genes <- function(cond1, cond2, df_deq = Dds_all, fdrlevel.de = 0.05){
  dds <- df_deq
  
  dds <- dds[, dds$condition == cond1 | 
               dds$condition == cond2 ]
  dds$condition <- droplevels(dds$condition)
  
  # Filter low expressed genes (at least 10 counts in total, or more than 0 in at least 50% of the samples)
  dds <- dds[rowSums(Counts(dds)) >= 10,]
  dds_good <- dds[rowSums(Counts(dds)==0) <= ncol(Counts(dds))*0.5,]
  
  dds_good <- estimateSizeFactors(dds_good)
  
  rowData(dds_good)$control_expr <- rowMeans(Counts(dds_good, normalize = TRUE)[,dds_good$condition == cond1])
  rowData(dds_good)$condition_expr <- rowMeans(Counts(dds_good, normalize = TRUE)[,dds_good$condition == cond2])
  
  # Perform diffential expression analysis
  dds_good <- DESeq(dds_good)
  res_good <- results(dds_good)
  
  res_good$control_expr <- rowData(dds_good)$control_expr
  res_good$condition_expr <- rowData(dds_good)$condition_expr
  
  dds_good <- dds_good[!is.na(res_good$padj)]
  res_good <- res_good %>% na.omit()
  
  res_good <- res_good[res_good$padj < fdrlevel.de,]
  res_good$status <- ifelse(res_good$control_expr < res_good$condition_expr,"upregulated","downregulated")
  res_good$comparison <- paste0(cond1,"_vs_",cond2)
  
  return(res_good)
  
}

# Perform the differential expresiion analysis
First_deseq <- Get_diff_genes(cond1 = "Multinucleated", cond2 = "Amoeba")

# Create a table with the information about the gene expression changes
First_deseq$status %>% table()

# Convert results to a dataframe
Results_MultinucleatedvsAmoeba <- as.data.frame(First_deseq)

# Save results to a text file
write.table(Results_MultinucleatedvsAmoeba, file = "DESeq_resultsCfragMvsA.txt", sep = "\t", row.names = TRUE)


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

# Create matrix with the genes and tpm of the different cell stages samples (Make sure your abundance files have the same gene order)
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



# Check distribution of difference of methylated gene body (fraction)
Cfrag_MultVsAmoeba$delta_6mAfraction_Multi_vs_Amoeba %>% hist( breaks = 30)

# Plot the histogram of the difference between the fraction of mApTs in genes
hist(Cfrag_MultVsAmoeba$delta_6mAfraction_Multi_vs_Amoeba, breaks = 30,
     main = "Creolimax      Multinucleated      vs      Amoeba
  Difference Between Genes mApT fraction",
     xlab = "Difference between mApT fraction",
     ylab = "Number of Genes")

# Filter the genes that have a difference of 0.015
Amoeba_hyper_ApTfrac <- Cfrag_MultVsAmoeba %>% filter(delta_6mAfraction_Multi_vs_Amoeba <= -0.015) %>% .$gene
Amoeba_hypo_ApTfrac <- Cfrag_MultVsAmoeba %>% filter(delta_6mAfraction_Multi_vs_Amoeba >= 0.015) %>% .$gene

# Check number of hyper and hypo methylated genes
length(Amoeba_hyper_ApTfrac)
length(Amoeba_hypo_ApTfrac) 

# Check the genes that share status in fraction and in total numbers
table(Amoeba_hyper_ApTsite %in% Amoeba_hyper_ApTfrac)
table(Amoeba_hypo_ApTsite %in% Amoeba_hypo_ApTfrac)


# Add a column to the matrix adding a name to the filtrated genes for clarity
Cfrag_tpm_matrix_by6mAfraction <- Cfrag_tpm_matrix %>% mutate(category = ifelse(gene_id %in% Amoeba_hyper_ApTfrac, "Amoeba",
                                                                                ifelse(gene_id %in% Amoeba_hypo_ApTfrac, "Multinucleated","Other"))) %>% 
  filter(category != "Other")

# Reshape matrix
Cfrag_tpm_matrix_by6mAfraction <- melt(Cfrag_tpm_matrix_by6mAfraction)

# Plot the genes selected by the difference of the fraction of mApTs against the methylome of each cell stage
Cfrag_tpm_matrix_by6mAfractiondiff_plot_gg <- ggplot(Cfrag_tpm_matrix_by6mAfraction, aes(x = category, y = log10(value), fill = variable)) +
  labs(x = "Hypermethylated genes", y = "log10(TPM)", title = "TPM vs 6mAT by methylated fraction", fill = "Transcriptome") +
  scale_fill_manual(values = c("Amoeba1" = "#FBE6C5FF", "Amoeba2" = "#F5BA98FF", "Amoeba3" = "#EE8A82FF", "Multinucleated1" = "#C8586CFF", "Multinucleated2" = "#9C3F5DFF", "Multinucleated3" = "#9C3F5DFF"), 
                    labels = c("Amoeba 1", "Amoeba 2", "Amoeba 3", "Multinucleated 1", "Multinucleated 2", "Multinucleated 3" )) +  
  geom_boxplot()


# ------------------------------------------------------------------------------
# Comparison of different expressed genes and different ApT methylated genes in cell stages
# ------------------------------------------------------------------------------

# Load DEseq results 
Deseq_MvsA <- fread("DESeq_resultsCfragMvsA.txt")

# Count the number of upregulated genes
Deseq_MvsA %>% filter(status == "upregulated") %>% nrow() 

# Count the number of upregulated genes
Deseq_MvsA%>% filter(status == "upregulated", log2FoldChange < -1.4  & padj < 0.01) %>% nrow()

# Count the number of downregulated genes
Deseq_MvsA%>% filter(status == "downregulated", log2FoldChange > 1.4  & padj < 0.01) %>% nrow()

# Get the IDs of upregulated genes
Deseq_MvsA_ids <- Deseq_MvsA %>% filter(status == "upregulated", log2FoldChange < -1.4  & padj < 0.01 ) %>% .$gene_id

# Get the IDs of downregulated genes
Down_deseq_MvsA_ids <- Deseq_MvsA %>% filter(status == "downregulated", log2FoldChange > 1.4  & padj < 0.01 ) %>% .$gene_id


# Create a cleaned dataframe
Cfrag_6m_df_clean <- Cfrag_MultVsAmoeba %>% dplyr::select(gene, methylated_ApT.x, mApT.x, methylated_ApT.y, mApT.y) %>%
  dplyr::rename(methylated_ApT_multinucleated = methylated_ApT.x, methylated_ApT_amoeba = methylated_ApT.y,
                mApT_multinucleated = mApT.x, mApT_amoeba = mApT.y)

# Add a new column with the category of each gene based on DESeq analysis
Cfrag_6m_df_clean <- Cfrag_6m_df_clean %>% mutate(category = ifelse(gene %in% Deseq_MvsA_ids, "Amoeba_up",
                                                                    ifelse(gene %in% Down_deseq_MvsA_ids, "Multinucleated_up","Other"))) %>%
  filter(category != "Other")



# ------------------------------------------------------------------------------
# Comparison between different species
# ------------------------------------------------------------------------------


################Pearson correlation value between sample/stages#################


# Load the values from each species and establish an order for plotting
Stages_treatment_species <- c("Ngru", "Acas", "Spun", "Cper", "Cfra")
Pearson_stages_values <- c(0.7532508, 0.9258102, 0.9459579, 0.8916315, 0.9367317) #Load the values obtained before in this script -> Cfra_AmoebavsMultinucleated_Pearson
Stages_species_order <- rev(c("Ngru", "Acas", "Spun", "Cper", "Cfra"))

# Generate a data frame with the values
ApT_stages_methylation_correlations <- data.frame(Stages_treatment_species, Pearson_stages_values)

# Transform to plot
df_long_ApT_Stages <- tidyr::gather(ApT_stages_methylation_correlations, key = "variable", value = "value", -Stages_treatment_species)

# Plot a bar chart of the values
Pearson_correlation_between_stages_plot_gg <- ggplot(df_long_ApT_Stages, aes(x = reorder(Stages_treatment_species, match(Stages_treatment_species, Stages_species_order)), y = value, fill = variable)) +
  geom_bar(stat = "identity", fill="#FFE099FF", position = position_dodge(width = 0.8)) +
  labs(title = "ApT Methylation per Sample Correlation", x= 'Species', y = "Pearson Correlation Value", fill = "Variable") +
  theme_light() +
  ylim(-0.002,1) +
  guides(fill = "none") +
  geom_hline(yintercept = 1, color = "grey65", linetype = "dashed", linewidth = 0.8)


################6mA gene body fraction difference across samples################


# Name your dataset (Each is build as explain before with Cfra as an example)
Spun_ZoosporesVSMultinucleated$Dataset <- "Spun"
Cfrag_MultVsAmoeba$Dataset <- "Cfra"
Cper_54hVS96H$Dataset <- "Cper"
Acas_CystVsAMoeba$Dataset <- "Acas"
Ngru_23Vs30$Dataset <- "Ngru"

# Select the difference in 6AmT fraction form each sample
Spun_selected_fraction <- Spun_ZoosporesVSMultinucleated %>%
  select(gene, Dataset, starts_with("delta_6mAf"))

Cfra_selected_fraction <- Cfrag_MultVsAmoeba %>%
  select(gene, Dataset, starts_with("delta_6mAf"))

Cper_selected_fraction <- Cper_54hVS96H %>%
  select(gene, Dataset, starts_with("delta_6mAf"))

Acas_selected_fraction <- Acas_CystVsAMoeba %>%
  select(gene, Dataset, starts_with("delta_6mAf"))

Ngru_selected_fraction <- Ngru_23Vs30 %>%
  select(gene, Dataset, starts_with("delta_6mAf"))

# Combine samples
Combined_delta_fraction<- bind_rows(Spun_selected_fraction, Cfra_selected_fraction, Cper_selected_fraction, Acas_selected_fraction, Ngru_selected_fraction)

# Arrange data for plotting
Long_data_fraction <- Combined_delta_fraction %>%
  pivot_longer(cols = starts_with("delta_"), names_to = "Delta_Type", values_to = "Delta_Value")

# Establish order of species for ploting
Long_data_fraction$Dataset <- factor(Long_data_fraction$Dataset, 
                                     levels = c("Ngru", "Acas", "Spun", "Cper", "Cfra"))

# Plot the 6mAT fraction difference across samples
mApTfrac_difference_between_samples_plot_gg <- ggplot(Long_data_fraction, aes(x = Delta_Value, y = Dataset, fill = Dataset)) +
  geom_density_ridges(scale = 1.3, alpha = 0.9) +
  xlim(-0.025, 0.025) +
  scale_fill_manual(values = c(
    "Cfra" = "#217A79FF",
    "Cper" = "#4C9B82FF",
    "Spun" = "#6CC08BFF",
    "Acas" = "#97E196FF",
    "Ngru" = "#D3F2A3FF"
  )) +
  theme_ridges() +
  labs(
    title = "6mAT Fraction Difference Between Samples",
    x = "6mAT Fraction Difference",
    y = "Number of Genes"
  ) +
  guides(fill = "none")


################Overlap of hyper and hypomethylated genes and highly transcribed genes across samples################


#Generate a data frame using highly transcribed and hypermethylated data from before 
Cfrag_Selected_HiglyTrans_HyperMeth <- Cfrag_6m_df_clean
Cfrag_Selected_HiglyTrans_HyperMeth$Hypermethylated_Sample <- ifelse(Cfra_Selected_HiglyTrans_HyperMeth$gene %in% amoeba_hyper_ApTfrac, 
                                                                    "Amoeba", 
                                                                    ifelse(Cfra_Selected_HiglyTrans_HyperMeth$gene %in% amoeba_hypo_ApTfrac, 
                                                                           "Multinucleated", 
                                                                           "None"))

#Get the number of overlaping genes
table(Cfra_Selected_HiglyTrans_HyperMeth$category, Cfra_Selected_HiglyTrans_HyperMeth$Hypermethylated_Sample)

# Build a dataframe with the number of overlaping genes of each species. Up-regulated and Downregulated are genes that not overlap with any hyper or hypo methylated gene.
df_overlaping_genes_per_frac <- data.frame( 
  Species = rep(c("Acas", "Ngru", "Spun", "Cfra", "Cper"), each = 6), 
  Condition = rep(c("Up-regulated", "Downregulated", "Up-Hypermethylated Overlap", "Up-Hypomethylated Overlap", "Down-Hypomethylated Overlap", "Down-Hypermethylated Overlap"), times = 5),  
  Count = c(1013, 605, 360, 81, 78, 99,  
            645, 736, 3, 0, 8, 5,
            1918, 1853, 60, 45, 86, 104,
            431, 52, 36, 12, 3, 8,
            849, 107, 138, 66, 12, 18)
)


#Get the percentages of each stage
df_overlaping_genes_per_frac <- df_overlaping_genes_per_frac %>%
  group_by(Species) %>%
  mutate(
    # Total for the Up-regulated group (Up-regulated + overlaps)
    Total_Upregulated = sum(Count[Condition %in% c("Up-regulated", 
                                              "Up-Hypermethylated Overlap", 
                                              "Up-Hypomethylated Overlap")]),
    # Total for the Downregulated group (Downregulated + overlaps)
    Total_Downregulated = sum(Count[Condition %in% c("Downregulated", 
                                              "Down-Hypomethylated Overlap", 
                                              "Down-Hypermethylated Overlap")]),
    # Compute Percentage relative to the group
    Percentage = case_when(
      Condition %in% c("Up-regulated", "Up-Hypermethylated Overlap", "Up-Hypomethylated Overlap") ~ 
        Count / Total_Upregulated * 100,  
      Condition %in% c("Downregulated", "Down-Hypomethylated Overlap", "Down-Hypermethylated Overlap") ~ 
        -Count / Total_Downregulated * 100
    )
  )


#Separate stages in positive and negative axis for plotting
df_overlaping_genes_per_frac$Adjusted_Count <- ifelse(df_overlaping_genes_per_frac$Condition %in% c("Up-regulated", "Up-Hypermethylated Overlap", "Up-Hypomethylated Overlap"), 
                                                      df_overlaping_genes_per_frac$Count, 
                                                      -df_overlaping_genes_per_frac$Count)

#Establish order for plotting
df_overlaping_genes_per_frac$Condition <- factor(df_overlaping_genes_per_frac$Condition, 
                                                 levels = c("Up-regulated", "Up-Hypermethylated Overlap", "Up-Hypomethylated Overlap", 
                                                            "Downregulated", "Down-Hypomethylated Overlap", "Down-Hypermethylated Overlap"))

#Plot as stacked bars the overlap between selected genes
Overlap_trans_meth_genes_stacked_plot_gg <- ggplot(df_overlaping_genes_per_frac, aes(x = reorder(Species, match(Species, Stages_species_order)), y = Adjusted_Count, fill = Condition)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(
    "Up-regulated" = "#B7E6A5FF",
    "Downregulated" = "#FFE099FF",
    "Up-Hypermethylated Overlap" = "#7CCBA2FF",
    "Up-Hypomethylated Overlap" = "#089099FF",
    "Down-Hypomethylated Overlap" = "#FF9C4CFF",
    "Down-Hypermethylated Overlap" = "#FF7200FF"
  )) +
  labs(title = "Overlap of Transcription and Methylated Genes",
       x = "Species",
       y = "Gene Number",
       fill = "Condition") +
  theme_minimal()



################################################################################
###H3K4me3 depletion vs 6mAT levels
################################################################################

# Prepare data frame with the 6mAT methylation percentages of each treatment group
Acas_Samples <- c("Amoeba", "Cyst", "DMSO 0.2%", "Piribedil 250uM", "DMSO 0.6%", "OICR-9429 300uM")
Global_6mAT_Perc <- c(3.47007, 3.16354, 3.53444, 3.40125, 2.56244, 2.47007 ) #This values are obtained as explain before from nanopore R9 sequencing files

# Establish the comparison groups between treatments
Comparison_Group <- c("Amoeba vs Cyst", "Amoeba vs Cyst", "Piribedil", "Piribedil", "OICR-9424", "OICR-9424")

# Create data frame
Acas_ApT_methylation_Samples <- data.frame(Acas_Samples, Global_6mAT_Perc, Comparison_Group)

# Plot a bar chart with the comparison of the 6mAT percentages per treatment 
Acas_Treatment_6mAT_Perc_Plot_gg <- ggplot(Acas_ApT_methylation_Samples, aes(x = Comparison_Group, y = Global_6mAT_Perc, fill = Acas_Samples)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) + 
  geom_text(aes(label = Acas_Samples), 
            position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("Amoeba" = "#ECDA9AFF", "Cyst" = "#EE4D5AFF", 
                               "DMSO 0.6%"= "#f5e6b3", "OICR-9429 300uM" = "#eb5461",
                               "DMSO 0.2%" = "#f7eabc", "Piribedil 250uM" = "#f5626f")) + # 
  labs(title = "Effect on H3k4m3 drugs on 6mAT%",
       x = "Comparison Groups",
       y = "Global 6mAT%") +
  theme_minimal()


########Unix command line code to plot 6mAT treatment levels vs H3K4me3 ########
#
## Use bigwigs obtained before from nanopore R9 sequencing files to plot against H3K4me3 narrow peaks bed file using Deeptools2
#
## Compute matrix using the Center as reference point
#computeMatrix reference-point --referencePoint center -out Acas_AmoebavsCyst.narrowPeak.R9mods.dtmatrix \
#-S Acas.R9_6mA.methyl.ApT.bigwig Acas_Cyst.R9_6mA.methyl.ApT.bigwig \
#-R Acas_H3K4me3.merge.dedup_peaks.narrowPeak.bed \
#-a 2000 -b 2000 -bs 10 -p 20
#
## Plot profile using the computed matrix
#plotProfile -m Acas_AmoebavsCyst.narrowPeak.R9mods.dtmatrix -o Acas_AmoebavsCyst.h3K4m3narrowPeak.R9mods_center.pdf \
#--color green blue --perGroup --samplesLabel Amoeba Cyst --yMax 26
#
## Repeat with the rest of the samples
#computeMatrix reference-point --referencePoint center -out Acas_OICRvsControl.narrowPeak.R9mods.dtmatrix \
#-S Acas_300OICR.R9_6mA.methyl.ApT.bigwig Acas_06DMSO.R9_6mA.methyl.ApT.bigwig \
#-R Acas_H3K4me3.merge.dedup_peaks.narrowPeak.bed \
#-a 2000 -b 2000 -bs 10 -p 20
#
#plotProfile -m Acas_OICRvsControl.narrowPeak.R9mods.dtmatrix -o Acas_OICRvsControl.h3K4m3narrowPeak.R9mods_center.pdf \
#--color orange red --perGroup --samplesLabel "OICR-9429 300uM" "DMSO 0.6%" --yMax 26
#
#computeMatrix reference-point --referencePoint center -out Acas_PiribedilvsControl.narrowPeak.R9mods.dtmatrix \
#-S Acas_250Piribedil.R9_6mA.methyl.ApT.bigwig Acas_02DMSO.R9_6mA.methyl.ApT.bigwig \
#-R Acas_H3K4me3.merge.dedup_peaks.narrowPeak.bed \
#-a 2000 -b 2000 -bs 10 -p 20
#
#plotProfile -m Acas_PiribedilvsControl.narrowPeak.R9mods.dtmatrix -o Acas_PiribedilvsControl.h3K4m3narrowPeak.R9mods_center.pdf \
#--color pink purple --perGroup --samplesLabel "Piribedil 250uM" "DMSO 0.2%" --yMax 26
#
