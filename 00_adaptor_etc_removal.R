# Adaptor, Primer, Pad, and Linker Removal

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
