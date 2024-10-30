library(dplyr)
library(seqR)
library(stringr)
library(kmerFilters)
library(bigstep)

sequence_dat <- readRDS("data/sequence_dat.RDS")

sequence_dat <- sequence_dat %>% 
  mutate(rep = str_sub(protein,-1,-1))

sequence_dat_rep1 <- sequence_dat %>% 
  filter(rep == "1") %>% 
  mutate(id = 1:n())

# kmer space

all_three_gaps <- expand.grid(0L:6, 0L:6)
gaps_shorter_than6 <- all_three_gaps[rowSums(all_three_gaps) <= 6, ]

k_vector_raw <- c(1, 2, 2, 2, 2, 2, rep(3, nrow(gaps_shorter_than6)))
kmer_gaps_raw <- c(list(NULL, 
                        NULL, c(1), c(2), c(3), c(4)),
                   lapply(1L:nrow(gaps_shorter_than6), function(i)
                     unlist(gaps_shorter_than6[i, ], use.names = FALSE)))

kmers <- count_multimers(
  sequence_dat_rep1[["sequence"]],
  k_vector = k_vector_raw,
  kmer_gaps = kmer_gaps_raw,
  with_kmer_counts = FALSE,
  batch_size = 4)

#####

# example train test split

set.seed(1)

test_ids <- sequence_dat_rep1 %>% 
  group_by(positive) %>% 
  sample_frac(0.3) %>% 
  pull(id)

kmers_test <- kmers[test_ids, ]
kmers_train <- kmers[-test_ids, ]

target <- sequence_dat_rep1 %>% 
  pull(positive) %>% 
  as.numeric()

target_train <- target[-test_ids]
target_test <- target[test_ids]


dat <- prepare_data(target_train, as.matrix(kmers_train))

dat <- reduce_matrix(dat, minpv = 0.3)

res_mbic2 <- stepwise(dat, crit = "mbic2")[["model"]]






