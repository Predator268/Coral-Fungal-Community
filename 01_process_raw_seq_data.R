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
