
library(dplyr)
library(seqinr)
library(biogram)
library(kmerFilters)


aa_sequences <- read.fasta("data/AMPBenchmark_public.fasta")

is_positive <- unlist(sapply(aa_sequences, function(ith_seq) {
  grepl("AMP=1", attr(ith_seq, "name"))
}))

sequences_dat <- data.frame(
  sequence = sapply(aa_sequences, function(x) paste0(x, collapse = ""),
                    USE.NAMES = FALSE),
  positive = unlist(sapply(aa_sequences, function(ith_seq) {
    grepl("AMP=1", attr(ith_seq, "name"))
  }))
)

protein <- rownames(sequences_dat)

sequences_dat <- cbind(protein = protein, sequences_dat)

rownames(sequences_dat) <- NULL

saveRDS(sequences_dat, "data/sequence_dat.RDS")



