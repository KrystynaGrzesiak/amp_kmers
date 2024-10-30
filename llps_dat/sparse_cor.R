
sparse_cor <- function(mbic_chosen, kmers_train, candidates, threshold) {
  
  n <- nrow(kmers_train)
  candidates_kmers <- kmers_train[, candidates]
  candidates_sums <- Matrix::colSums(candidates_kmers)
  
  for(y_name in mbic_chosen) {
    
    xy <- Matrix::colSums(candidates_kmers * candidates_kmers[, y_name])
    
    nom <- n * xy - candidates_sums * candidates_sums[y_name]
    denom <- sqrt(n * candidates_sums - candidates_sums^2) *
      sqrt(n * candidates_sums[y_name] - candidates_sums[y_name]^2)
    
    correlations <- nom / denom
    
    chosen_kmers <- c(chosen_kmers, names(which(correlations >= threshold)))
  }
  chosen_kmers
}
