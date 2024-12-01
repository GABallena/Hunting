if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(GenomicRanges)
library(AnnotationHub)
library(GenomicFeatures)

# Load the genome
chr22_seq <- Hsapiens$chr22

# Question 1: GC Content of chr22 in hg19
gc_content_chr22 <- function(sequence) {
  gc <- sum(letterFrequency(sequence, letters = c("G", "C")))
  total <- sum(letterFrequency(sequence, letters = c("A", "T", "G", "C")))
  gc / total
}
gc_chr22 <- gc_content_chr22(chr22_seq)
cat("Question 1: GC Content of chr22:", gc_chr22, "\n")

# Question 2: Mean GC Content in NarrowPeak Regions
narrowPeaks <- import("path_to_narrowPeak_file.bed", format = "BED")
gc_content_regions <- function(start, end) {
  region_seq <- subseq(chr22_seq, start, end)
  gc_content_chr22(region_seq)
}
gc_values <- mapply(gc_content_regions, start(narrowPeaks), end(narrowPeaks))
mean_gc_content <- mean(gc_values)
cat("Question 2: Mean GC Content in NarrowPeak regions:", mean_gc_content, "\n")

# Question 3: Correlation Between GC Content and SignalValue
signal_values <- mcols(narrowPeaks)$score
cor_gc_signal <- cor(gc_values, signal_values)
cat("Question 3: Correlation between GC content and signalValue:", cor_gc_signal, "\n")

# Question 4: Correlation Between SignalValue and Average fc.signal
hub <- AnnotationHub()
fc_signal_data <- query(hub, c("H3K27me3", "fold-change", "H1"))
fc_signal <- hub[["AH12345"]]  # Replace with actual data ID from AnnotationHub

mean_fc_signal <- mapply(function(start, end) {
  mean(fc_signal[start:end])
}, start(narrowPeaks), end(narrowPeaks))
cor_signal_fc <- cor(signal_values, mean_fc_signal)
cat("Question 4: Correlation between signalValue and average fc.signal:", cor_signal_fc, "\n")

# Question 5: Number of Bases with fc.signal >= 1
high_signal_bases <- sum(fc_signal[seqnames(fc_signal) == "chr22"] >= 1)
cat("Question 5: Bases on chr22 with fc.signal >= 1:", high_signal_bases, "\n")

# Question 6: Regions with E003 Signal <= 0.5 and E055 Signal >= 2
e003_signal <- fc_signal  # Replace with actual E003 signal data
e055_signal <- fc_signal  # Replace with actual E055 signal data

diff_regions <- which(e003_signal <= 0.5 & e055_signal >= 2)
cat("Question 6: Regions where E003 <= 0.5 and E055 >= 2:", length(diff_regions), "\n")

# Question 7: CpG Islands - Observed-to-Expected Ratio
observed_to_expected_cpg <- function(sequence) {
  cpg_obs <- countPattern("CG", sequence)
  c_freq <- letterFrequency(sequence, "C") / width(sequence)
  g_freq <- letterFrequency(sequence, "G") / width(sequence)
  cpg_exp <- c_freq * g_freq * width(sequence)
  cpg_obs / cpg_exp
}
cpg_islands <- import("path_to_CpG_islands.bed", format = "BED")
o_to_e_ratios <- sapply(seq_along(cpg_islands), function(i) {
  observed_to_expected_cpg(subseq(chr22_seq, start(cpg_islands[i]), end(cpg_islands[i])))
})
average_o_to_e_ratio <- mean(o_to_e_ratios)
cat("Question 7: Average observed-to-expected ratio of CpG islands:", average_o_to_e_ratio, "\n")

# Question 8: Count of TATA Boxes on chr22
tata_box_count <- countPattern("TATAAA", chr22_seq) + countPattern("TTTATA", chr22_seq)
cat("Question 8: Number of TATA boxes on chr22:", tata_box_count, "\n")

# Question 9: TATA Boxes in Promoters of Transcripts
txdb <- makeTxDbFromUCSC(genome = "hg19", tablename = "knownGene")
promoters <- promoters(txdb, upstream = 900, downstream = 100)
promoters_chr22 <- subset(promoters, seqnames(promoters) == "chr22")
tata_promoters <- sum(mapply(function(start, end, strand) {
  promoter_seq <- subseq(chr22_seq, start, end)
  countPattern("TATAAA", promoter_seq) > 0
}, start(promoters_chr22), end(promoters_chr22), strand(promoters_chr22)))
cat("Question 9: Promoters with TATA box:", tata_promoters, "\n")

# Question 10: Overlapping Promoters
overlapping_bases <- sum(width(reduce(promoters_chr22, ignore.strand = TRUE)) > 1)
cat("Question 10: Bases part of more than one promoter:", overlapping_bases, "\n")
