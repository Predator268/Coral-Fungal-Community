############################################
#           00_adaptor_removal.R           #
############################################

# Adaptor and Primer Removal

# Raw data is used as input to the 00_adaptor_etc_removal.R script to remove adaptors and primers from sequences.
# This outputs N_filtered folder, Adaptor_Removed folder and Adaptor_Primer_Removed folder, all inside the raw data folder.
# The folders are fed to the script in the above order and the final folder is the one used for analysis.

library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")


# Declaring path to raw data
path <- "../Coral_fastq/Pa"
list.files(path)

fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))


# Declaring adaptor and primmers
FWD_adaptor <- "AATGATACGGCGACCACCGAGATCTACAC"  
REV_adaptor <- "AATGATACGGCGACCACCGAGATCTACAC"

FWD_primer <- "CTTGGTCATTTAGAGGAAGTAA"  
REV_primer <- "GCTGCGTTCTTCATCGATGC"

# Getting all orientations of each sequence
allOrients <- function(sequence) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(sequence)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD_adaptor.orients <- allOrients(FWD_adaptor)
REV_adaptor.orients <- allOrients(REV_adaptor)
FWD_adaptor.orients

FWD_primer.orients <- allOrients(FWD_primer)
REV_primer.orients <- allOrients(REV_primer)
FWD_primer.orients


# Pre-filtering to remove ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)


### Removal Process ###

#Using cutadapt for trimming
cutadapt <- "C:/Users/User/AppData/Roaming/Python/Python38/Scripts/cutadapt"
system2(cutadapt, args = "--version")

## Adaptor Removal ##
# Counting presence of adaptors before removal
sequenceHits <- function(sequence, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(sequence, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD_adaptor.ForwardReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnFs.filtN[[10]]), 
      FWD_adaptor.ReverseReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnRs.filtN[[10]]), 
      REV_adaptor.ForwardReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnFs.filtN[[10]]), 
      REV_adaptor.ReverseReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnRs.filtN[[10]]))

# Removing adaptors
path.cut_adaptor <- file.path(path, "Adaptor_Removed")
if(!dir.exists(path.cut_adaptor)) dir.create(path.cut_adaptor)
fnFs.cut_adaptor <- file.path(path.cut_adaptor, basename(fnFs))
fnRs.cut_adaptor <- file.path(path.cut_adaptor, basename(fnRs))

FWD_adaptor.RC <- dada2:::rc(FWD_adaptor)
REV_adaptor.RC <- dada2:::rc(REV_adaptor)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_adaptor, "-a", REV_adaptor.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_adaptor, "-A", FWD_adaptor.RC) 

# Run Cutadapt
outputStatsAdaptor <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "-o", fnFs.cut_adaptor[i], "-p", fnRs.cut_adaptor[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
  }
)
cat(outputStatsAdaptor, file="./Output/cutadapt_adaptor_trimming_stats.txt", sep="\n", append = FALSE)

# Counting adaptors after removal
rbind(FWD_adaptor.ForwardReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnFs.cut_adaptor[[10]]), 
      FWD_adaptor.ReverseReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnRs.cut_adaptor[[10]]), 
      REV_adaptor.ForwardReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnFs.cut_adaptor[[10]]), 
      REV_adaptor.ReverseReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnRs.cut_adaptor[[10]]))

## Primer Removal
#Counting Primers before removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.cut_adaptor[[10]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.cut_adaptor[[10]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.cut_adaptor[[10]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.cut_adaptor[[10]]))

# Removing primers
path.cut_primer <- file.path(path, "Primer_Adaptor_Removed")
if(!dir.exists(path.cut_primer)) dir.create(path.cut_primer)
fnFs.cut_primer <- file.path(path.cut_primer, basename(fnFs))
fnRs.cut_primer <- file.path(path.cut_primer, basename(fnRs))

FWD_primer.RC <- dada2:::rc(FWD_primer)
REV_primer.RC <- dada2:::rc(REV_primer)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_primer, "-a", REV_primer.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_primer, "-A", FWD_primer.RC) 

# Run Cutadapt
outputStatsPrimer <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "-o", fnFs.cut_primer[i], "-p", fnRs.cut_primer[i], # output files
                               fnFs.cut_adaptor[i], fnRs.cut_adaptor[i])) # input files
  }
)
cat(outputStatsPrimer, file="./Output/cutadapt_primer_trimming_stats.txt", sep="\n", append = FALSE)

# Counting primers after removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.cut_primer[[10]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.cut_primer[[10]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.cut_primer[[10]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.cut_primer[[10]]))



# Forward and reverse fastq filenames have the format:
cutFs_raw <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_raw <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
cut_raw <- sort(list.files(path, pattern = "fastq.gz", full.names = TRUE))

path.filtN <- file.path(path, "N_filtered")
cutFs_filtN <- sort(list.files(path.filtN, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_filtN <- sort(list.files(path.filtN, pattern = "_2.fastq.gz", full.names = TRUE))
cut_filtN <- sort(list.files(path.filtN, pattern = "fastq.gz", full.names = TRUE))

cutFs_adaptor <- sort(list.files(path.cut_adaptor, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_adaptor <- sort(list.files(path.cut_adaptor, pattern = "_2.fastq.gz", full.names = TRUE))
cut_adaptor <- sort(list.files(path.cut_adaptor, pattern = "fastq.gz", full.names = TRUE))

cutFs_primer <- sort(list.files(path.cut_primer, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_primer <- sort(list.files(path.cut_primer, pattern = "_2.fastq.gz", full.names = TRUE))
cut_primer <- sort(list.files(path.cut_primer, pattern = "fastq.gz", full.names = TRUE))


# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs_primer, get.sample.name))
head(sample.names)


## Quality Profiles ##
# Raw
plotQualityProfile(cutFs_raw[6:11])
plotQualityProfile(cutRs_raw[6:11])

# Pre-filtered
plotQualityProfile(cutFs_filtN[6:11])
plotQualityProfile(cutRs_filtN[6:11])

# Adaptor Removed
plotQualityProfile(cutFs_adaptor[6:11])
plotQualityProfile(cutRs_adaptor[6:11])

# Primer Removed
plotQualityProfile(cutFs_primer[6:11])
plotQualityProfile(cutRs_primer[6:11])


## Comparing files and sample properties for each step
library(seqTools)
packageVersion("seqTools")
library(microseq)
packageVersion("microseq")

samples <- sort(list.files(path, pattern = "fastq.gz", full.names = FALSE))


# Comparing number of reads for each sample during each step

raw_seqs <- numeric(length(cut_raw))
prefiltered_seqs <- numeric(length(cut_raw))
adaptorremoved_seqs <- numeric(length(cut_raw))
primerremoved_seqs <- numeric(length(cut_raw))
fractionOfRawDataInFinal_seqs <- primerremoved_seqs / raw_seqs

for (j in 1:length(samples)){
  raw_seqs[j] <- nrow(readFastq(cut_raw[j]))
  prefiltered_seqs[j] <- nrow(readFastq(cut_filtN[j]))
  adaptorremoved_seqs[j] <- nrow(readFastq(cut_adaptor[j]))
  primerremoved_seqs[j] <- nrow(readFastq(cut_primer[j]))
}
fractionOfRawDataInFinal_seqs <- primerremoved_seqs / raw_seqs


sample_summaries <- data.frame(samples, raw_seqs, prefiltered_seqs, adaptorremoved_seqs, primerremoved_seqs, fractionOfRawDataInFinal_seqs, stringsAsFactors = FALSE)
write.csv(sample_summaries, file = "Number of Sequence Summary.csv")
summary(sample_summaries)
hist(sample_summaries$raw_seqs)
hist(sample_summaries$primerremoved_seqs)










############################################
#        01_process_raw_seq_data.R         #
############################################

# Process raw sequences into phyloseq object for analyses

# Load dada2 and prep ####
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")


# File parsing - For this, we will use only the forward illumina reads
path <- "../Coral_fastq/Pa/Primer_Adaptor_Removed" # CHANGE to the directory containing your demultiplexed fastq files
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present
fns <- sort(list.files(file.path(path), full.names = TRUE, pattern = "_1.fastq.gz"))

sample.names <- unlist(map(strsplit(basename(fns), "_"), 1)) 
sample.names

# visualize a couple of fwd read quality profiles to help select reasonable filtration parameters
plotQualityProfile(fns[20:24])


# Filter and trim ####
filts <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))

out <- filterAndTrim(fns, filts, # fnRs, filtRs,
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=4) # On Windows set multithread=FALSE
saveRDS(out,"./Output/Pa/out.RDS")
#out <- readRDS("./Output/Pa/out.RDS")

# sanity check  comparison of before and after filtration
plotQualityProfile(c(fns[1:2],filts[1:2]))

# LEARN ERROR RATES ####
# Since some samples may have had zero reads pass QC, reassign filts
filts <- sort(list.files(filtpath, full.names = TRUE))
errF <- learnErrors(filts, multithread=TRUE, MAX_CONSIST = 20)
saveRDS(errF,"./Output/Pa/errF.RDS")
#errF <- readRDS("./Output/Pa/errF.RDS")

# sanity check for error model
plotErrors(errF, nominalQ=TRUE)

# DEREPLICATION ####
derep <- derepFastq(filts, verbose=TRUE)


# Name the derep-class objects by the sample names
# If some samples were removed (no reads passed QC), reassign sample.names
if(length(derep) != length(sample.names)){
  sample.names <- unlist(map(strsplit(basename(filts), "_filt"), 1))
}
names(derep) <- sample.names

# SAMPLE INFERRENCE ####
dadaFs <- dada(derep, err=errF, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, pool = "pseudo")
saveRDS(dadaFs,"./Output/Pa/dadaFs.RDS")
dadaFs <- readRDS("./Output/Pa/dadaFs.RDS")

# Make a sequence table ####
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


# REMOVE CHIMERAS ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# reassign "out" to remove any missing reads
out = out[as.data.frame(out)$reads.out > 0,]

# TRACK READS THROUGH PIPELINE ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out[,1], sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "nonchim")
rownames(track) <- sample.names
track = as.data.frame(track)
track$total.loss.proportion = (track[,1]-track[,2])/track[,1]
head(track)
write.csv(track, file = "./Output/Pa/read_counts_at_each_step.csv", row.names = TRUE)


# Save intermediate seqtable object
saveRDS(seqtab.nochim, "./Output/Pa/seqtab.nochim.RDS")
# Resuming Point
#seqtab.nochim <- readRDS("./Output/Pa/seqtab.nochim.RDS")

# Removing objects not required ahead to conserve RAM because some of the next steps are memory intensive
rm(dadaFs)
rm(derep)
rm(errF)
rm(out)
rm(seqtab)

# IMPORT METADATA ####
meta = as.data.frame(read_delim("./Pa_Meta.csv",delim = ","))
row.names(meta) <- as.character(meta$SampleID)
meta = meta[order(meta$SampleID),]
if(!identical(row.names(meta),row.names(seqtab.nochim))){
  keepers <- row.names(meta) %in% row.names(seqtab.nochim)
  meta <- meta[keepers,]  
}


# Remove all seqs with fewer than 100 nucleotides ####
df.track <- as.data.frame(Sample <- row.names(meta))
for (i in 1:nrow(meta)){
  df.track$beforeLengthFilter[i] <- sum(seqtab.nochim[i,])
}
keeper_esvs <- nchar(names(as.data.frame(seqtab.nochim))) > 99
seqtab.nochim <- seqtab.nochim[,keeper_esvs]

# ASSIGN TAXONOMY ####
taxa <- assignTaxonomy(seqtab.nochim, "./sh_general_release_dynamic_s_10.10.2017.fasta", multithread=4,verbose=TRUE)

# Save intermediate files
saveRDS(seqtab.nochim, file = "./Output/Pa/dada2_seqtable.RDS")
saveRDS(taxa, file = "./Output/Pa/RDP_Taxonomy_from_dada2.RDS")

# re-load point
seqtab.nochim <- readRDS("./Output/Pa/dada2_seqtable.RDS")
taxa <- readRDS("./Output/Pa/RDP_Taxonomy_from_dada2.RDS")


# Hand off to Phyloseq ####
otu <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
tax <- tax_table(taxa)
met <- sample_data(meta)
row.names(met) <- row.names(meta)


ps <- phyloseq(otu,met,tax)


# Find controlsamples (extraction negatives) and clean ####
contams <- decontam::isContaminant(otu_table(ps),neg = ps@sam_data$Location == "Blank",method = "prevalence")

# Remove contaminants 
for (i in 1:nrow(meta)){
  df.track$beforeDecontam[i] <- sum(ps@otu_table[i,])
}
ps.noncontam <- prune_taxa(!contams$contaminant, ps)
for (i in 1:nrow(meta)){
  df.track$afterDecontam[i] <- sum(ps.noncontam@otu_table[i,])
}

# REMOVE NON-FUNGI and empty samples/taxa ####
ps.noncontam <- subset_taxa(ps.noncontam, Kingdom == "k__Fungi")
for (i in 1:nrow(meta)){
  df.track$fungiOnly[i] <- sum(ps.noncontam@otu_table[i,])
}
ps.noncontam <- subset_taxa(ps.noncontam, taxa_sums(ps.noncontam) > 0)
ps.noncontam <- subset_samples(ps.noncontam, sample_sums(ps.noncontam) > 0)
for (i in 1:nrow(meta)){
  df.track$nonEmpty[i] <- sum(ps.noncontam@otu_table[i,])
}

# Save DNA sequences apart from rownames (from subsetted ps object)
seqs <- taxa_names(ps)
seqs <- DNAStringSet(seqs)
saveRDS(seqs,"./Output/Pa/ASV_reference_sequences.RDS")


pretty_names <- paste("FungalASV",1:length(taxa_names(ps)),":",
                      tax_table(ps)[,2],
                      tax_table(ps)[,3],
                      tax_table(ps)[,4],
                      tax_table(ps)[,5],
                      tax_table(ps)[,6],
                      tax_table(ps)[,7], sep="_") %>%
  str_remove("k__") %>% str_remove("p__") %>% str_remove("c__") %>% str_remove("o__") %>% str_remove("f__") %>% str_remove("g__") %>% str_remove("s__") %>%
  str_replace(pattern = "_:_",replacement = ": ")

df <- data.frame(TaxaName=pretty_names,Sequence=taxa_names(ps))
saveRDS(df,"./Output/Pa/SequenceNames_and_Taxonomy.RDS")

ps.noncontam <- prune_samples(ps.noncontam@sam_data$Location != "Blank", ps.noncontam)

# Save RDS object for Phyloseq
saveRDS(ps.noncontam, file = "./Output/Pa/clean_phyloseq_object.RDS")

# Final Filtering
# load processed data
ps = readRDS(file = "./Output/Pa/clean_phyloseq_object.RDS")

# quick peek at ps object
sample_names(ps)
rank_names(ps)
sample_variables(ps)
otu_table(ps)[1:5, 1:5]
tax_table(ps)[1:5, 1:6]

# Change sequences to unique IDs to make viewing easier
seqs_ITS = taxa_names(ps)
names(seqs_ITS) <- 1:length(seqs_ITS) # save seqs and IDs combination
taxa_names(ps) <- names(seqs_ITS)
tax_table(ps)[1:5, 1:6]


## Removing less prevalent sequences
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))

prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps, "Phylum"))

#  Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps)

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps1 = prune_taxa(keepTaxa, ps)

nBlanks <- 6
startPoint <- 1+nBlanks
df.track$afterPrevFilter <- numeric(nrow(df.track))
for (j in startPoint:nrow(meta)){ 
  df.track$afterPrevFilter[j] <- sum(ps1@otu_table[(j-nBlanks),])
}

# Remove samples not associated with an island
ps1 = prune_samples(ps1@sam_data$Island != "", ps1)
df.track$afterLocationFilter_Final <- numeric(nrow(df.track))
for (j in startPoint:nrow(meta)){ 
  df.track$afterLocationFilter_Final[j] <- sum(ps1@otu_table[(j-nBlanks),])
}

# Final ps object
saveRDS(ps1, file = "./Output/Pa/final_phyloseq_object.RDS")

write.csv(df.track, file = "./Output/Pa/FinalTrackReads.csv", row.names = TRUE)










############################################
#              02_analysis.R               #
############################################

# Loading Packages
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(broom)
library(RColorBrewer)
library(RgoogleMaps)
library(metagMisc)
library(shiny)
library(viridis)
library(svglite)
library(gdtools)
library(VennDiagram)
library(indicspecies)
library(geosphere)
library(ecodist)

# Color palette
pal = c("#6b5456","#ec8d1b","#6abf2a","#8b53b7","#70acbe","#01c95b","#c00014","#31332f","#f7d000","#abba00")
col_island <- viridis(9)
names(col_island) <- c("Hantu", "Jong", "Kusu", "Raffles Lighthouse", "Semakau", "Sisters", "St. John", "Sultan Shoal", "TPT")

theme_set(theme_bw())

#### Analysis ####
ps1 <- readRDS("./Output/Pa/final_phyloseq_object.RDS")
taxtable.update <- as.data.frame(ps1@tax_table)
for(j in 1:ncol(taxtable.update)){
  taxtable.update[,j] <- as.character(taxtable.update[,j])
  for(i in 1:nrow(taxtable.update)){
    if(!is.na(taxtable.update[i, j])){
      textLabel <- strsplit(as.character(taxtable.update[i, j]), split = '_')[[1]][3:length(strsplit(as.character(taxtable.update[i, j]), split = '_')[[1]])]
      taxtable.update[i, j] <- paste(textLabel, collapse='_')
    }
  }
}


# Data Prep for FUNGuild (Actual run done through python script Guilds_v1.1.py)
taxafinal <- taxtable.update
taxafinal$taxonomy <- taxafinal$Species
for(i in 1:nrow(taxafinal)){
  taxafinal$taxonomy[i] <- paste(taxafinal[i,-ncol(taxafinal)], collapse = ' ')
}
OTUTable <- as.data.frame(ps1@otu_table)
OTUTable <- as.data.frame(t(OTUTable))
OTUTable$OTU_ID <- row.names(OTUTable)
OTUTable$taxonomy <- taxafinal$taxonomy
OTUTable <- OTUTable[,98:99]
#write.csv(OTUTable, "./Output/Pa/FUNGuild_OTU_Table.csv")


# Normalize (relative abundance) ####
ps1ra <- transform_sample_counts(ps1, function(otu){otu/sum(otu)})



### Rarefaction Curves ###
grp <- factor(ps1@sam_data$Island)
cols <- col_island[grp]
rarefaction_curve <- rarecurve(ps1@otu_table, step = 20, col = cols, label = FALSE)
Nmax <- sapply(rarefaction_curve, function(x) max(attr(x, "Subsample")))
Smax <- sapply(rarefaction_curve, max)
svg(filename = "./Output/Pa/Analysis/Rarefaction Curves.svg", )
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Number of Sequences",
     ylab = "Number of ASVs", type = "n",
     main = "Rarefaction Curves")
for (i in seq_along(rarefaction_curve)) {
  N <- attr(rarefaction_curve[[i]], "Subsample")
  lines(N, rarefaction_curve[[i]], col = cols[i])
}
legend(9000, 7, legend = names(col_island), col = col_island, lty = 1, cex = 0.8, box.lty = 1)
dev.off()





### Relative Abundance Bar Plots by Location and Structure ###

# Creating a dataframe of the otu table which also includes the Island variables
combineddf <- as.data.frame(ps1@otu_table)
combineddf$Island <- ps1@sam_data$Island
##############
# Creating a phyloseq object where samples are grouped and merged according to Island
island_otu <- combineddf %>%
  group_by(Island) %>%
  summarise_all(.funs=sum)
island_otu <- as.data.frame(island_otu)
row.names(island_otu) <- c("Samples_Hantu", "Samples_Jong", "Samples_Kusu", "Samples_Raffles", "Samples_Semakau", "Samples_Sisters", "Samples_John", "Samples_Sultan", "Samples_TPT")
island_otu <- island_otu[,-1]
island_meta <- data.frame(Island = c("Hantu", "Jong", "Kusu", "Raffles Lighthouse", "Semakau", "Sisters", "St. John", "Sultan Shoal", "TPT"))
row.names(island_meta) <- NULL
row.names(island_meta) <- c("Samples_Hantu", "Samples_Jong", "Samples_Kusu", "Samples_Raffles", "Samples_Semakau", "Samples_Sisters", "Samples_John", "Samples_Sultan", "Samples_TPT")
speciesList <- taxtable.update$Species
genusList <- taxtable.update$Genus
genusSpeciesList <- speciesList
for(i in 1:length(speciesList)){
  if(!is.na(speciesList[i])){
    genusSpeciesList[i] <- paste(genusList[i], speciesList[i], sep = ' ')
  }
}
taxtable.update$Genus_species <- genusSpeciesList

island_ps <- phyloseq(otu_table(island_otu, taxa_are_rows=FALSE), 
                      sample_data(island_meta), 
                      tax_table(ps1@tax_table))

# Calculating Relative abundances and plotting bar plots according to Location and Structure
ps_Island_ra <- transform_sample_counts(island_ps, function(otu){otu/sum(otu)})


# Defining function to make bar charts without black lines separating samples. Based on phyloseq function "plot_bar".
simple_plot_bar = function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, 
                            facet_wrap = NULL) {
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack")
  p = p + theme_bw() + theme(axis.text=element_text(size=15), axis.title=element_text(size=17,face="bold"), 
                             axis.text.x = element_text(angle = 70, hjust = 1))
  p = p + labs(y = "Relative Abundance")
  p = p + guides(guide_legend(ncol = 1), fill = guide_legend(ncol = 3))
  if (!is.null(facet_wrap)) {
    p <- p + facet_wrap(facet_wrap, nrow = 1) + theme(strip.text = element_text(size=15))
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

# Making stacked bar charts for relative abundance of taxa

# According to Phylum
ra_Phylum_barplot_island <- simple_plot_bar(ps_Island_ra, x="Island", fill="Phylum") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Phylum)))
ggsave(ra_Phylum_barplot_island, filename = "./Output/Pa/Analysis/Relative Abundance of Ohylum by Island - Bar Plot.svg", dpi=300, width = 12, height = 10)

ra_Phylum_barplot_all <- simple_plot_bar(ps1ra, x="SampleID", fill="Phylum") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Phylum)))
ggsave(ra_Phylum_barplot_all, filename = "./Output/Pa/Analysis/Relative Abundance of Phylum by Sample - Bar Plot.svg", dpi=300, width = 25, height = 10)

ra_Phylum_barplot_combined <- simple_plot_bar(ps1ra, x="sample_Species", fill="Phylum") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Phylum)))
ggsave(ra_Phylum_barplot_combined, filename = "./Output/Pa/Analysis/Relative Abundance of Phylum Combined - Bar Plot.svg", dpi=300, width = 12, height = 10)

# According to Class
ra_Class_barplot_island <- simple_plot_bar(ps_Island_ra, x="Island", fill="Class") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Class)))
ggsave(ra_Class_barplot_island, filename = "./Output/Pa/Analysis/Relative Abundance of Class by Island - Bar Plot.svg", dpi=300, width = 12, height = 10)

ra_Class_barplot_all <- simple_plot_bar(ps1ra, x="SampleID", fill="Class") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Class)))
ggsave(ra_Class_barplot_all, filename = "./Output/Pa/Analysis/Relative Abundance of Class by Sample - Bar Plot.svg", dpi=300, width = 25, height = 10)

ra_Class_barplot_combined <- simple_plot_bar(ps1ra, x="sample_Species", fill="Class") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Class)))
ggsave(ra_Class_barplot_combined, filename = "./Output/Pa/Analysis/Relative Abundance of Class Combined - Bar Plot.svg", dpi=300, width = 12, height = 10)

# According to Order
ra_Order_barplot_island <- simple_plot_bar(ps_Island_ra, x="Island", fill="Order") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Order)))
ggsave(ra_Order_barplot_island, filename = "./Output/Pa/Analysis/Relative Abundance of Order by Island - Bar Plot.svg", dpi=300, width = 12, height = 10)

ra_Order_barplot_all <- simple_plot_bar(ps1ra, x="SampleID", fill="Order") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Order)))
ggsave(ra_Order_barplot_all, filename = "./Output/Pa/Analysis/Relative Abundance of Order by Sample - Bar Plot.svg", dpi=300, width = 25, height = 10)

ra_Order_barplot_combined <- simple_plot_bar(ps1ra, x="sample_Species", fill="Order") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Order)))
ggsave(ra_Order_barplot_combined, filename = "./Output/Pa/Analysis/Relative Abundance of Order Combined - Bar Plot.svg", dpi=300, width = 12, height = 10)

# According to Family
ra_Family_barplot_island <- simple_plot_bar(ps_Island_ra, x="Island", fill="Family") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Family)))
ggsave(ra_Family_barplot_island, filename = "./Output/Pa/Analysis/Relative Abundance of Family by Island - Bar Plot.svg", dpi=300, width = 12, height = 10)

ra_Family_barplot_all <- simple_plot_bar(ps1ra, x="SampleID", fill="Family") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Family)))
ggsave(ra_Family_barplot_all, filename = "./Output/Pa/Analysis/Relative Abundance of Family by Sample - Bar Plot.svg", dpi=300, width = 25, height = 10)

ra_Family_barplot_combined <- simple_plot_bar(ps1ra, x="sample_Species", fill="Family") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Family)))
ggsave(ra_Family_barplot_combined, filename = "./Output/Pa/Analysis/Relative Abundance of Family Combined - Bar Plot.svg", dpi=300, width = 12, height = 10)

# According to Genus
ra_Genus_barplot_island <- simple_plot_bar(ps_Island_ra, x="Island", fill="Genus") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Genus)))
ggsave(ra_Genus_barplot_island, filename = "./Output/Pa/Analysis/Relative Abundance of Genus by Island - Bar Plot.svg", dpi=300, width = 12, height = 10)

ra_Genus_barplot_all <- simple_plot_bar(ps1ra, x="SampleID", fill="Genus") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Genus)))
ggsave(ra_Genus_barplot_all, filename = "./Output/Pa/Analysis/Relative Abundance of Genus by Sample - Bar Plot.svg", dpi=300, width = 25, height = 10)

ra_Genus_barplot_combined <- simple_plot_bar(ps1ra, x="sample_Species", fill="Genus") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = sort(unique(taxtable.update$Genus)))
ggsave(ra_Genus_barplot_combined, filename = "./Output/Pa/Analysis/Relative Abundance of Genus Combined - Bar Plot.svg", dpi=300, width = 12, height = 10)

# According to Species
sorter <- c(sort(unique(taxtable.update$Species), index.return = TRUE)$ix + 1, 1)

ra_Species_barplot_island <- simple_plot_bar(ps_Island_ra, x="Island", fill="Species") + theme(legend.text=element_text(size=rel(1.3)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = unique(taxtable.update$Genus_species)[sorter])
ggsave(ra_Species_barplot_island, filename = "./Output/Pa/Analysis/Relative Abundance of Species by Island - Bar Plot.svg", dpi=300, width = 12, height = 10)

ra_Species_barplot_all <- simple_plot_bar(ps1ra, x="SampleID", fill="Species") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = unique(taxtable.update$Genus_species)[sorter])
ggsave(ra_Species_barplot_all, filename = "./Output/Pa/Analysis/Relative Abundance of Species by Sample - Bar Plot.svg", dpi=300, width = 25, height = 10)

ra_Species_barplot_combined <- simple_plot_bar(ps1ra, x="sample_Species", fill="Species") + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80)) + scale_fill_discrete(labels = unique(taxtable.update$Genus_species)[sorter])
ggsave(ra_Species_barplot_combined, filename = "./Output/Pa/Analysis/Relative Abundance of Species Combined - Bar Plot.svg", dpi=300, width = 12, height = 10)




### Heatmap ###
melted_ps <- psmelt(ps1ra)


## By Phylum
factor_melted_ps_by_sample <- unique(melted_ps$SampleID) #factoring by sample
dataframe <- data.frame()  #creating dataframe to fill in information

for (i in factor_melted_ps_by_sample) {                                    #For each sample,
  sub_ps <- melted_ps[melted_ps$SampleID == i,]  #subset dataframe.
  
  per_phylum_abundance <- aggregate(Abundance ~ Phylum, data = sub_ps, sum)  #Aggregate abundance based on each phylum
  
  per_phylum_abundance$SampleID <- sub_ps$SampleID[1:nrow(per_phylum_abundance)] #Adding sample name to each phylum row
  per_phylum_abundance$Island <- sub_ps$Island[1:nrow(per_phylum_abundance)]
  
  dataframe <- rbind(dataframe, per_phylum_abundance)                       #store this in a dataframe for each row
}

#Sorting dataframe in order of island
ordered_df <- dataframe[order(dataframe$Island),]

#Plotting
table(ordered_df$Island) / sum(unique(ordered_df$Phylum) == unique(ordered_df$Phylum))
heatplot <- ggplot(ordered_df, aes(reorder(SampleID, -desc(Island)), reorder(Phylum, desc(Phylum)))) +
  geom_tile(aes(fill = Abundance)) + 
  labs(y = "Phylum", x = "Samples") + 
  theme(axis.text.x = element_blank(), axis.title = element_text(size = 20), axis.text.y = element_text(size = 13), legend.text = element_text(size = 13), legend.title = element_text(size = 15)) + 
  scale_fill_gradient(low = "#FFFFFF", high = "#680000") +
  geom_vline(xintercept = c(58, 117, 177), alpha = 0.2)

#Adding labels
y.min <- 2.5
y.max <- 3.5
y.mid <- (y.min + y.max) / 2
heatplot <- heatplot +
  annotate("rect", xmin = 0.5, xmax = 5.5, ymin = y.min, ymax = y.max,
           alpha = 0.9, fill = col_island[1]) +
  annotate("rect", xmin = 5.5, xmax = 14.5, ymin = y.min, ymax = y.max,
           alpha = 0.95, fill = col_island[2]) +
  annotate("rect", xmin = 14.5, xmax = 34.5, ymin = y.min, ymax = y.max,
           alpha = 0.95, fill = col_island[3]) +
  annotate("rect", xmin = 34.5, xmax = 54.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[4]) + 
  annotate("rect", xmin = 54.5, xmax = 65.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[5]) + 
  annotate("rect", xmin = 65.5, xmax = 85.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[6]) + 
  annotate("rect", xmin = 85.5, xmax = 86.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[7]) + 
  annotate("rect", xmin = 86.5, xmax = 87.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[8]) + 
  annotate("rect", xmin = 87.5, xmax = 97.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[9]) + 
  annotate("text", x = 3, y = y.mid, label = names(col_island[1]), size = 7, srt = 90) +
  annotate("text", x = 10, y = y.mid, label = names(col_island[2]), size = 7, srt = 90) +
  annotate("text", x = 24.5, y = y.mid, label = names(col_island[3]), size = 7, srt = 90) +
  annotate("text", x = 44.5, y = y.mid, label = names(col_island[4]), size = 7, srt = 90) +
  annotate("text", x = 60, y = y.mid, label = names(col_island[5]), size = 7, srt = 90) +
  annotate("text", x = 75.5, y = y.mid, label = names(col_island[6]), size = 7, srt = 90) +
  annotate("text", x = 86, y = y.mid, label = names(col_island[7]), size = 4, srt = 90) +
  annotate("text", x = 87, y = y.mid, label = names(col_island[8]), size = 4, srt = 90) +
  annotate("text", x = 92.5, y = y.mid, label = names(col_island[9]), size = 7, srt = 90)
heatplot
ggsave(heatplot, filename = "./Output/Pa/Analysis/Heatmap of Phylum Grouped by Island.svg", dpi = 300, width = 18, height = 10)


## By Class
factor_melted_ps_by_sample <- unique(melted_ps$Sample) #factoring by sample
dataframe <- data.frame()  #creating dataframe to fill in information

for (i in factor_melted_ps_by_sample) {                                    #For each sample,
  sub_ps <- melted_ps[melted_ps$SampleID == i,]  #subset dataframe.
  
  per_class_abundance <- aggregate(Abundance ~ Class, data = sub_ps, sum)  #Aggregate abundance based on each class
  
  per_class_abundance$SampleID <- sub_ps$SampleID[1:nrow(per_class_abundance)] #Adding sample name to each class row
  per_class_abundance$Island <- sub_ps$Island[1:nrow(per_class_abundance)]
  
  dataframe <- rbind(dataframe, per_class_abundance)                       #store this in a dataframe for each row
}

#Sorting dataframe in order of structure, sub-ordered by location
ordered_df <- dataframe[order(dataframe$Island),]

#Plotting
table(ordered_df$Island) / sum(unique(ordered_df$Class) == unique(ordered_df$Class))
heatplot <- ggplot(ordered_df, aes(reorder(SampleID, -desc(Island)), reorder(Class, desc(Class)))) +
  geom_tile(aes(fill = Abundance)) + 
  labs(y = "Class", x = "Samples") + 
  theme(axis.text.x = element_blank(), axis.title = element_text(size = 20), axis.text.y = element_text(size = 10), legend.text = element_text(size = 9), legend.title = element_text(size = 10)) + 
  scale_fill_gradient(low = "#FFFFFF", high = "#680000") +
  geom_vline(xintercept = c(58, 117, 177), alpha = 0.2)

#Adding labels
y.min <- 4.5
y.max <- 6
y.mid <- (y.min + y.max) / 2
heatplot <- heatplot +
  annotate("rect", xmin = 0.5, xmax = 5.5, ymin = y.min, ymax = y.max,
           alpha = 0.9, fill = col_island[1]) +
  annotate("rect", xmin = 5.5, xmax = 14.5, ymin = y.min, ymax = y.max,
           alpha = 0.95, fill = col_island[2]) +
  annotate("rect", xmin = 14.5, xmax = 34.5, ymin = y.min, ymax = y.max,
           alpha = 0.95, fill = col_island[3]) +
  annotate("rect", xmin = 34.5, xmax = 54.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[4]) + 
  annotate("rect", xmin = 54.5, xmax = 65.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[5]) + 
  annotate("rect", xmin = 65.5, xmax = 85.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[6]) + 
  annotate("rect", xmin = 85.5, xmax = 86.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[7]) + 
  annotate("rect", xmin = 86.5, xmax = 87.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[8]) + 
  annotate("rect", xmin = 87.5, xmax = 97.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[9]) + 
  annotate("text", x = 3, y = y.mid, label = names(col_island[1]), size = 7, srt = 90) +
  annotate("text", x = 10, y = y.mid, label = names(col_island[2]), size = 7, srt = 90) +
  annotate("text", x = 24.5, y = y.mid, label = names(col_island[3]), size = 7, srt = 90) +
  annotate("text", x = 44.5, y = y.mid, label = names(col_island[4]), size = 7, srt = 90) +
  annotate("text", x = 60, y = y.mid, label = names(col_island[5]), size = 7, srt = 90) +
  annotate("text", x = 75.5, y = y.mid, label = names(col_island[6]), size = 7, srt = 90) +
  annotate("text", x = 86, y = y.mid, label = names(col_island[7]), size = 4, srt = 90) +
  annotate("text", x = 87, y = y.mid, label = names(col_island[8]), size = 4, srt = 90) +
  annotate("text", x = 92.5, y = y.mid, label = names(col_island[9]), size = 7, srt = 90)
heatplot
ggsave(heatplot, filename = "./Output/Pa/Analysis/Heatmap of Class Grouped by Island.svg", dpi = 300, width = 18, height = 10)


## By Species
factor_melted_ps_by_sample <- unique(melted_ps$Sample) #factoring by sample
dataframe <- data.frame()  #creating dataframe to fill in information

for (i in factor_melted_ps_by_sample) {                                    #For each sample,
  sub_ps <- melted_ps[melted_ps$SampleID == i,]  #subset dataframe.
  
  per_species_abundance <- aggregate(Abundance ~ Species, data = sub_ps, sum)  #Aggregate abundance based on each species
  
  per_species_abundance$SampleID <- sub_ps$SampleID[1:nrow(per_species_abundance)] #Adding sample name to each species row
  per_species_abundance$Island <- sub_ps$Island[1:nrow(per_species_abundance)]
  
  dataframe <- rbind(dataframe, per_species_abundance)                       #store this in a dataframe for each row
}

#Sorting dataframe in order of structure, sub-ordered by location
ordered_df <- dataframe[order(dataframe$Island),]

#Plotting
table(ordered_df$Island) / sum(unique(ordered_df$Species) == unique(ordered_df$Species))
genus.species.labels <- unique(taxtable.update$Genus_species)[rev(sorter)]
genus.species.labels <- genus.species.labels[-1]
heatplot <- ggplot(ordered_df, aes(reorder(SampleID, -desc(Island)), reorder(Species, desc(Species)))) +
  geom_tile(aes(fill = Abundance)) + 
  labs(y = "Species", x = "Samples") + 
  theme(axis.text.x = element_blank(), axis.title = element_text(size = 20), axis.text.y = element_text(size = 10), legend.text = element_text(size = 9), legend.title = element_text(size = 10)) + 
  scale_fill_gradient(low = "#FFFFFF", high = "#680000") +
  geom_vline(xintercept = c(58, 117, 177), alpha = 0.2) +
  scale_y_discrete(labels = genus.species.labels)

#Adding labels
y.min <- 21.5
y.max <- 29
y.mid <- (y.min + y.max) / 2
heatplot <- heatplot +
  annotate("rect", xmin = 0.5, xmax = 5.5, ymin = y.min, ymax = y.max,
           alpha = 0.9, fill = col_island[1]) +
  annotate("rect", xmin = 5.5, xmax = 14.5, ymin = y.min, ymax = y.max,
           alpha = 0.95, fill = col_island[2]) +
  annotate("rect", xmin = 14.5, xmax = 34.5, ymin = y.min, ymax = y.max,
           alpha = 0.95, fill = col_island[3]) +
  annotate("rect", xmin = 34.5, xmax = 54.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[4]) + 
  annotate("rect", xmin = 54.5, xmax = 65.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[5]) + 
  annotate("rect", xmin = 65.5, xmax = 85.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[6]) + 
  annotate("rect", xmin = 85.5, xmax = 86.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[7]) + 
  annotate("rect", xmin = 86.5, xmax = 87.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[8]) + 
  annotate("rect", xmin = 87.5, xmax = 97.5, ymin = y.min, ymax = y.max,
           alpha = 1, fill = col_island[9]) + 
  annotate("text", x = 3, y = y.mid, label = names(col_island[1]), size = 7, srt = 90) +
  annotate("text", x = 10, y = y.mid, label = names(col_island[2]), size = 7, srt = 90) +
  annotate("text", x = 24.5, y = y.mid, label = names(col_island[3]), size = 7, srt = 90) +
  annotate("text", x = 44.5, y = y.mid, label = names(col_island[4]), size = 7, srt = 90) +
  annotate("text", x = 60, y = y.mid, label = names(col_island[5]), size = 7, srt = 90) +
  annotate("text", x = 75.5, y = y.mid, label = names(col_island[6]), size = 7, srt = 90) +
  annotate("text", x = 86, y = y.mid, label = names(col_island[7]), size = 4, srt = 90) +
  annotate("text", x = 87, y = y.mid, label = names(col_island[8]), size = 4, srt = 90) +
  annotate("text", x = 92.5, y = y.mid, label = names(col_island[9]), size = 7, srt = 90)
heatplot
ggsave(heatplot, filename = "./Output/Pa/Analysis/Heatmap of Species Grouped by Island.svg", dpi = 300, width = 18, height = 10)




### Shannon Diversity Plots ###
div <- data.frame(Island = ps1ra@sam_data$Island,
                  Shannon = diversity(otu_table(ps1ra)))
Richness = colSums(decostand(otu_table(ps1ra), method = "pa"))
write.csv(div, file = "./Output/Pa/Analysis/Diversity_table.csv", quote = FALSE)
# By Island
div %>% group_by(Island) %>%
  summarise(N = n(), Mean = mean(Shannon))
ggplot(div, aes(x=Island, y=Shannon, fill = Island)) + 
  geom_boxplot() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 25), axis.title.y = element_text(vjust = 2), axis.title = element_text(size = 30), title = element_text(size = 17), legend.text = element_text(size = 25), legend.title = element_text(size = 30)) + 
  labs(y="Shannon Diversity") +
  scale_fill_manual(values = col_island)
ggsave(filename = "./Output/Pa/Analysis/Shannon_Diversity_by_Location.svg", dpi = 300, width = 12, height = 10)




### Ordination(s) ###
ps.nmds <- subset_taxa(ps1ra, !is.na(Genus)) # Getting rid off non identified taxa at genus level
ps.nmds <- prune_samples(sample_sums(ps.nmds) > 0, ps.nmds)

NMDS = ordinate(ps.nmds, method = "NMDS", distance = "bray", trymax = 100)
PCoA = ordinate(ps.nmds, method = "PCoA", distance = "bray", trymax = 100)

# Stress Plot
svg("./Output/Pa/Analysis/Full_NMDS_Stress_Plot.svg")
p_stress_full <- stressplot(NMDS, title("Stress Plot for NMDS"))
dev.off()

# NMDS with/without ellipses
NMDS2 = data.frame(NMDS1 = NMDS$points[,1], NMDS2 = NMDS$points[,2],group=ps.nmds@sam_data$Island)
NMDS2.mean=aggregate(NMDS2[,1:2],list(group=ps.nmds@sam_data$Island),mean)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame()
for(g in levels(NMDS2$group)){
  if(g != "St. John" & g != "Sultan Shoal")
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS2[NMDS2$group==g,],
                                                     veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                  ,group=g))
}
p_NMDS1 <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = group), size = 4) +
  ggtitle(paste("NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + theme_bw() + scale_color_manual(values = col_island) + 
  theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
        legend.text = element_text(size = 20)) + 
  guides(shape = guide_legend(override.aes = list(size = 5)), color = guide_legend(override.aes = list(size = 5))) + 
  labs(color = "Island")
p_NMDS1_ell <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = group), size = 4) +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=group), size=2, linetype=2) +
  ggtitle(paste("NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + theme_bw() + scale_color_manual(values = col_island) + 
  theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
        legend.text = element_text(size = 20)) + 
  guides(shape = guide_legend(override.aes = list(size = 5)), color = guide_legend(override.aes = list(size = 5))) + 
  labs(color = "Structure")
p_NMDS1_ell

# PCoA
p_PCoA <- plot_ordination(ps.nmds, PCoA, color = "Island") +  
  theme_bw() + scale_color_manual(values = col_island) + geom_point(size = 3) +
  theme(axis.title = element_text(size = 10), title = element_text(size = 12), 
        legend.text = element_text(size = 10)) + 
  guides(color = guide_legend(override.aes = list(size = 4)))

ggsave(p_NMDS1, filename = "./Output/Pa/Analysis/Full_NMDS_w_Island_colored.svg", dpi = 300, width = 12, height = 10)
ggsave(p_NMDS1_ell, filename = "./Output/Pa/Analysis/Full_NMDS_w_Island_colored_and_ellipses.svg", dpi = 300, width = 12, height = 10)
ggsave(p_PCoA, filename = "./Output/Pa/Analysis/Full_PCoA_w_Island_colored.svg", dpi = 300, width = 12, height = 10)




### PERMANOVA Test ###
ps.permanova <- subset_taxa(ps.nmds, !is.na(Genus))

otu = as.data.frame(otu_table(ps.permanova))
meta = as.data.frame(sample_data(ps.permanova))
df = data.frame(SampleID = meta$SampleID, Island = meta$Island)
# Island
permanova_Island <- adonis(otu ~ Island, data = df)

sink("./Output/Pa/Analysis/adonis_Island_table.txt")
noquote(print("PermANOVA Table:"))
permanova_Island
sink(NULL)

pairwise.adonis <- function(x,factors, sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  co = combn(unique(factors),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(co[1,elem],co[2,elem]),] ~ factors[factors %in% c(co[1,elem],co[2,elem])] , method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}
padonis_Island <- pairwise.adonis(otu,as.character(meta$Island))
sink("./Output/Pa/Analysis/adonis_Island_table.txt", append = TRUE)
noquote(print("Pairwise adonis between islands (Bonferroni corrected Pvalues):"))
padonis_Island
sink(NULL)




### Mantel Test ###
ps.mantel <- ps1ra

## Extracting Longitude and Latitude data
meta <- as.data.frame(ps.mantel@sam_data)
meta$lon <- sapply(strsplit(meta$GPS, " "), `[`, 2)
meta$lon <- as.double(sapply(strsplit(meta$lon, "E"), `[`, 1), length=10)
meta$lat <- as.double(sapply(strsplit(meta$GPS, "N"), `[`, 1), length=10)

geo_full <- data.frame(meta$lon, meta$lat)

## Preparing asv tables
otu_full <- as.data.frame(ps.mantel@otu_table)

## Making distance matrices
# Adundance data frames - bray curtis dissimilarity
dist.otu_full <- vegdist(otu_full, method = "bray")

# Geographic data frame - haversie distance
d.geo_full <- distm(geo_full, fun = distHaversine)

dist.geo_full <- as.dist(d.geo_full)

## Running Mantel Test
# Abundance vs Geographic
abund_geo_full <- vegan::mantel(dist.otu_full, dist.geo_full, method = "spearman", permutations = 9999, na.rm = TRUE)

## Saving Output
sink("./Output/Pa/Analysis/Mantel_Test.txt")
noquote("Mantel Test on All Samples")
abund_geo_full
sink(NULL)




### Multiple Regression on distance matrices ###
dist_MRM <- MRM(dist.otu_full ~ dist.geo_full,  nperm = 9999)
dist_MRM

sink("./Output/Pa/Analysis/MRM_Table.txt")
print("Bray-Curtis distance regressed against spatial distance (Multiple regression on matrices) (All Samples):")
print(dist_MRM)
sink(NULL)

#################################