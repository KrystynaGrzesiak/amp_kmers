
library(biogram)
library(seqR)
library(ranger)
library(bigstep)
library(boot)
library(pROC)
library(ggplot2)
library(caret)
library(dplyr)
library(glmnet)
library(tidyr)
library(SLOPE)
library(ggcorrplot)

inv_logit <- function(xb) exp(xb)/(1 + exp(xb))

df <- read.table("llps_dat/llps_ml_dataset.csv", sep = "\t", header = TRUE)
train_ids <- df[["Fold"]] == "Train"
test_ids <- df[["Fold"]] == "Test"


### kmer space

all_three_gaps <- expand.grid(0L:6, 0L:6)
gaps_shorter_than6 <- all_three_gaps[rowSums(all_three_gaps) <= 6, ]

k_vector_raw <- c(1, 2, 2, 2, 2, 2, rep(3, nrow(gaps_shorter_than6)))
kmer_gaps_raw <- c(list(NULL, 
                        NULL, c(1), c(2), c(3), c(4)),
                   lapply(1L:nrow(gaps_shorter_than6), function(i)
                     unlist(gaps_shorter_than6[i, ], use.names = FALSE)))

kmers <- count_multimers(
  df[["Full.seq"]],
  k_vector = k_vector_raw,
  kmer_gaps = kmer_gaps_raw,
  with_kmer_counts = FALSE,
  batch_size = 4)

######

y <- !(df[["Datasets"]] %in% c("PDB", "DisProt"))

train_y <- as.numeric(y[train_ids])
train_x <- data.frame((as.matrix(kmers[train_ids, ])))

test_y <- as.numeric(y[test_ids])
test_x <- data.frame((as.matrix(kmers[test_ids, ])))



###### mbic2 normal and logistic

dat <- prepare_data(train_y, train_x)

dat_reduced <- dat %>%
  reduce_matrix(minpv = 0.2)

res_mbic2 <- dat_reduced %>%
  stepwise(crit = "mbic2")


saveRDS(res_mbic2[["model"]], "llps_dat/mbic_vars.RDS")


# To się liczy bardzo długo i nie daje innego wyniku
# dat <- prepare_data(train_y, train_x, type = "logistic")
# 
# dat_reduced_logistic <- dat %>% 
#   reduce_matrix(minpv = 0.2)
# 
# res_mbic2_logistic <- dat_reduced %>% 
#   stepwise(crit = "mbic2")
# 
# saveRDS(res_mbic2_logistic, "mbic2_results_logistic.RDS")

# res_mbic2 <- readRDS("llps_dat/mbic2_results.RDS")


chosen_kmers <- readRDS("llps_dat/mbic_vars.RDS")
train_filtered_mbic2 <- train_x[, chosen_kmers]

test_y <- as.numeric(y[test_ids])
test_x <- data.frame((as.matrix(kmers[test_ids, chosen_kmers])))

cost_function <- function(r, pi = 0) { mean(abs(r - pi)) }

# model bez interakcji

logistic_model <- glm(train_y~., data = train_filtered_mbic2, family = "binomial")
predicted <- predict(logistic_model, test_x, type = "response")
rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)

cv_result <- cv.glm(data = cbind(train_y, train_filtered_mbic2), 
                    glmfit = logistic_model, K = 10)
cv_result$delta 

confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))

# model z interakcjami

logistic_model <- glm(train_y~(.)^2, data = train_filtered_mbic2, family = "binomial")
predicted <- predict(logistic_model, test_x, type = "response")
rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)
cv_result <- cv.glm(data = cbind(train_y, train_filtered_mbic2), 
                    glmfit = logistic_model, K = 10)
cv_result$delta 

confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))


########################## pełne dane

dat <- prepare_data(train_y, train_x, type = "logistic")

dat_reduced <- dat %>%
  reduce_matrix(minpv = 0.1) # to się długo liczy

# saveRDS(dat_reduced$candidates, "llps_dat/reduced_02.RDS")

candidates <- readRDS("llps_dat/reduced_02.RDS")

train_x_reduced <- train_x[, candidates]
test_x_reduced <- data.frame((as.matrix(kmers[test_ids, candidates])))

# res_lasso <- glmnet(train_x_reduced, train_y, family = "binomial") # to się liczy chwile

saveRDS(res_lasso, "llps_dat/lasso_results.RDS")
# cv_lasso_res <- cv.glmnet(train_x_reduced, train_y, family = "binomial") # za długie to jest

inv_logit <- function(xb) exp(xb)/(1 + exp(xb))

predicted <- predict.glmnet(res_lasso, newx = as.matrix(test_x_reduced), type = "response")
predicted <- apply(predicted, 2, inv_logit)

dt_coefs <- data.frame(as.matrix(coef.glmnet(res_lasso))) 

lambda_df <- data.frame(path = names(res_lasso$a0), lambda = res_lasso$lambda)

dt_coefs <- dt_coefs %>% 
  mutate(beta = rownames(dt_coefs)) %>% 
  tidyr::gather(path, estimator, -beta) %>% 
  left_join(lambda_df, by = "path")

dt_coefs %>% 
  filter(beta %in% chosen_kmers) %>% 
  ggplot(aes(x = lambda, y = estimator, group_by(beta), col = beta)) +
  geom_line()


predicted_single_lambda <- predicted[, 40]
rocobj <- roc(test_y, predicted_single_lambda)
ggroc(rocobj)
auc(test_y, predicted_single_lambda)

sum(coef.glmnet(res_lasso)[, 50] != 0)

aucs <- sapply(1:ncol(predicted), function(ith_lambda) {
  predicted_single_lambda <- predicted[, ith_lambda]
  auc(test_y, predicted_single_lambda)
})

data.frame(auc = aucs,
           lambda = res_lasso$lambda) %>% 
  ggplot(aes(x = lambda, y = auc)) +
  geom_line() +
  scale_x_reverse()


confusionMatrix(as.factor(as.numeric(predicted_single_lambda > 0.5)), as.factor(test_y))

##### SLOPE

res_slope <- SLOPE(train_x_reduced, train_y, family = "binomial", lambda = "bh", alpha = 0.001)

saveRDS(res_slope, "llps_dat/res_slope.RDS")

predicted <- predict(res_slope, as.matrix(test_x_reduced), type = "response")

rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)

################################################################################
################################################################################
#################################### 17-10-2024 ################################
################################################################################
################################################################################

# 2. LASSO, sLOPE, ridge, elastic net
# - na zbiorze z interakcjami 
# - na zbiorze z korelacjami z mbic + ewentualnie interakcje

# dat <- prepare_data(train_y, train_x)
# 
# dat_reduced <- dat %>%
#   reduce_matrix(minpv = 0.2)
# 
# res_mbic2 <- dat_reduced %>%
#   stepwise(crit = "mbic2")

get_interactions_design <- function(y, x) {
  lm_model <- lm(y~ (.)^2, data = data.frame(x))
  model.matrix(lm_model)
}

source("llps_dat/sparse_cor.R")

chosen_kmers <- readRDS("llps_dat/mbic_vars.RDS")

train_x_reduced <- as.matrix(train_x[, chosen_kmers])
test_x_reduced <- as.matrix(test_x[, chosen_kmers])

design_matrix_train <- get_interactions_design(train_y, train_x_reduced)
design_matrix_test <- get_interactions_design(test_y, test_x_reduced)

ggcorrplot(cor(design_matrix_train))

### z interakcjami 

# glm

logistic_model <- glm(train_y~., data = as.data.frame(design_matrix_train), family = "binomial")
predicted <- predict(logistic_model, as.data.frame(design_matrix_test), type = "response")
rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)

# LASSO

lasso_cv <- cv.glmnet(design_matrix_train, train_y, family = "binomial")

lambda_min <- lasso_cv$lambda.min

lasso_res <- glmnet(design_matrix_train, train_y, family = "binomial",
                    lambda = lambda_min, thresh = 10^(-10), maxit = 10^(8))

predicted <- predict.glmnet(lasso_res, newx = as.matrix(design_matrix_test), type = "response")
predicted <- as.vector(inv_logit(predicted))

rocobj <- roc(test_y, as.vector(predicted))
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))



# elastic net

elastic_net <- glmnet(design_matrix_train, train_y, family = "binomial", 
                      alpha = 0.5, lambda = lambda_min, thresh = 10^(-10), 
                      maxit = 10^(8))

predicted <- predict.glmnet(elastic_net, newx = as.matrix(design_matrix_test), type = "response")
predicted <- as.vector(apply(predicted, 2, inv_logit))

rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))


# ridge

ridge <- glmnet(design_matrix_train, train_y, family = "binomial", 
                lambda = lambda_min, alpha = 0, thresh = 10^(-10), maxit = 10^(8))

predicted <- predict.glmnet(ridge, newx = as.matrix(design_matrix_test), type = "response")
predicted <- as.vector(apply(predicted, 2, inv_logit))

rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))

## slope

res_slope <- SLOPE(design_matrix_train, train_y, family = "binomial", 
                   lambda = "bh", alpha = 0.01)


predicted <- predict(res_slope, as.matrix(design_matrix_test), type = "response")

rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))


### z korelacjami

# #### To się długo liczyło bardzio
# correlated_kmers <- sparse_cor(mbic_chosen = chosen_kmers, 
#                            train_x, 
#                            dat[["candidates"]],
#                            threshold = 0.8)
# 
# saveRDS(unique(correlated_kmers), "llps_dat/correlated_kmers.RDS")


chosen_kmers <- readRDS("llps_dat/mbic_vars.RDS")

# kmery skorelowane na poziomie 0.8
correlated_kmers <- readRDS("llps_dat/correlated_kmers.RDS")

# to jest to samo

###################################

# 3. Dużo większa liczba zmiennych, np. 1000, 3000, 5000 (wybrane testami brzegowymi)
# - zapuścić wszystko tak jak było


# dat <- prepare_data(train_y, train_x)
# 
# dat_reduced <- dat %>%
#   reduce_matrix(minpv = 1)
# 
# saveRDS(dat_reduced$candidates, "llps_dat/ranking.RDS")

ranking <- readRDS("llps_dat/ranking.RDS")


# 1000

candidates <- ranking[1:1000]

train_x_reduced <- as.matrix(train_x[, candidates])
test_x_reduced <- as.matrix(test_x[, candidates])

# lassso

lasso_cv <- cv.glmnet(train_x_reduced, train_y, family = "binomial")

lambda_min <- lasso_cv$lambda.min

lasso_res <- glmnet(train_x_reduced, train_y, family = "binomial",
                    lambda = lambda_min, thresh = 10^(-10), maxit = 10^(8))

predicted <- predict.glmnet(lasso_res, newx = as.matrix(test_x_reduced), type = "response")
predicted <- as.vector(inv_logit(predicted))

rocobj <- roc(test_y, as.vector(predicted))
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))


# elastic net

elastic_net <- glmnet(train_x_reduced, train_y, family = "binomial", 
                      alpha = 0.5, lambda = lambda_min, thresh = 10^(-10), 
                      maxit = 10^(8))

predicted <- predict.glmnet(elastic_net, newx = as.matrix(test_x_reduced), type = "response")
predicted <- as.vector(apply(predicted, 2, inv_logit))

rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))


# ridge

ridge <- glmnet(train_x_reduced, train_y, family = "binomial", 
                lambda = lambda_min, alpha = 0, thresh = 10^(-10), maxit = 10^(8))

predicted <- predict.glmnet(ridge, newx = as.matrix(test_x_reduced), type = "response")
predicted <- as.vector(apply(predicted, 2, inv_logit))

rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))

## slope

res_slope <- SLOPE(train_x_reduced, train_y, family = "binomial", 
                   lambda = "bh", alpha = 0.006)


predicted <- predict(res_slope, as.matrix(test_x_reduced), type = "response")

rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))

################

# 3000

candidates <- ranking[1:3000]

train_x_reduced <- as.matrix(train_x[, candidates])
test_x_reduced <- as.matrix(test_x[, candidates])

# lassso

lasso_cv <- cv.glmnet(train_x_reduced, train_y, family = "binomial")

lambda_min <- lasso_cv$lambda.min

lasso_res <- glmnet(train_x_reduced, train_y, family = "binomial",
                    lambda = lambda_min, thresh = 10^(-10), maxit = 10^(8))

predicted <- predict.glmnet(lasso_res, newx = as.matrix(test_x_reduced), type = "response")
predicted <- as.vector(inv_logit(predicted))

rocobj <- roc(test_y, as.vector(predicted))
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))


# elastic net

elastic_net <- glmnet(train_x_reduced, train_y, family = "binomial", 
                      alpha = 0.5, lambda = lambda_min, thresh = 10^(-10), 
                      maxit = 10^(8))

predicted <- predict.glmnet(elastic_net, newx = as.matrix(test_x_reduced), type = "response")
predicted <- as.vector(apply(predicted, 2, inv_logit))

rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))


# ridge

ridge <- glmnet(train_x_reduced, train_y, family = "binomial", 
                lambda = lambda_min, alpha = 0, thresh = 10^(-10), maxit = 10^(8))

predicted <- predict.glmnet(ridge, newx = as.matrix(test_x_reduced), type = "response")
predicted <- as.vector(apply(predicted, 2, inv_logit))

rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))

## slope

res_slope <- SLOPE(train_x_reduced, train_y, family = "binomial", 
                   lambda = "bh", alpha = 0.006)


predicted <- predict(res_slope, as.matrix(test_x_reduced), type = "response")

rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))




################

# 5000

candidates <- ranking[1:5000]

train_x_reduced <- as.matrix(train_x[, candidates])
test_x_reduced <- as.matrix(test_x[, candidates])

# lassso

lasso_cv <- cv.glmnet(train_x_reduced, train_y, family = "binomial")

lambda_min <- lasso_cv$lambda.min

lasso_res <- glmnet(train_x_reduced, train_y, family = "binomial",
                    lambda = lambda_min, thresh = 10^(-10), maxit = 10^(8))

predicted <- predict.glmnet(lasso_res, newx = as.matrix(test_x_reduced), type = "response")
predicted <- as.vector(inv_logit(predicted))

rocobj <- roc(test_y, as.vector(predicted))
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))


# elastic net

elastic_net <- glmnet(train_x_reduced, train_y, family = "binomial", 
                      alpha = 0.5, lambda = lambda_min, thresh = 10^(-10), 
                      maxit = 10^(8))

predicted <- predict.glmnet(elastic_net, newx = as.matrix(test_x_reduced), type = "response")
predicted <- as.vector(apply(predicted, 2, inv_logit))

rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))


# ridge

# ridge <- glmnet(train_x_reduced, train_y, family = "binomial", 
#                 lambda = lambda_min, alpha = 0, thresh = 10^(-10), maxit = 10^(8))
# 
# predicted <- predict.glmnet(ridge, newx = as.matrix(test_x_reduced), type = "response")
# predicted <- as.vector(apply(predicted, 2, inv_logit))
# 
# rocobj <- roc(test_y, predicted)
# ggroc(rocobj)
# auc(test_y, predicted)
# confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))

## slope

res_slope <- SLOPE(train_x_reduced, train_y, family = "binomial", 
                   lambda = "bh", alpha = 0.006)


predicted <- predict(res_slope, as.matrix(test_x_reduced), type = "response")

rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))




candidates <- ranking[1:10000]

train_x_reduced <- as.matrix(train_x[, candidates])
test_x_reduced <- as.matrix(test_x[, candidates])


res_slope <- SLOPE(train_x_reduced, train_y, family = "binomial", 
                   lambda = "bh", alpha = 0.006)


predicted <- predict(res_slope, as.matrix(test_x_reduced), type = "response")

rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))




candidates <- ranking[1:100]

train_x_reduced <- as.matrix(train_x[, candidates])
test_x_reduced <- as.matrix(test_x[, candidates])


res_slope <- SLOPE(train_x_reduced, train_y, family = "binomial", 
                   lambda = "bh", alpha = 0.006)



predicted <- predict(res_slope, as.matrix(test_x_reduced), type = "response")

rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))


#######################################################################

library(xgboost)

#### xgboost

chosen_kmers <- readRDS("../mbic_vars.RDS")
train_x_reduced <- as.matrix(train_x[, chosen_kmers])
test_x_reduced <- as.matrix(test_x[, chosen_kmers])

res_xgboost <- xgboost(data = train_x_reduced, label = train_y, max.depth = 2, eta = 1, 
        nthread = 2, nrounds = 2, objective = "binary:logistic")

predicted <- predict(res_xgboost, test_x_reduced)

rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))


#### random forest

model_rf_small <- ranger(train_y ~ .,
                         data = train_x_reduced,
                         probability = TRUE)

predicted <- predict(model_rf_small, test_x_reduced)[["predictions"]][, 1]
rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))

#### superlearner

library(SuperLearner)

sl <- SuperLearner(Y = train_y, X = train_x_reduced, family = binomial(),
                  SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"))


predicted <- as.vector(unlist(predict(sl, test_x_reduced)$pred))


################################################################################
################################################################################
#################################### 24-10-2024 ################################
################################################################################
################################################################################


get_interactions_design_3 <- function(y, x) {
  lm_model <- lm(y~ (.)^3, data = data.frame(x))
  model.matrix(lm_model)
}


chosen_kmers <- readRDS("llps_dat/mbic_vars.RDS")

train_x_reduced <- as.matrix(train_x[, chosen_kmers])
test_x_reduced <- as.matrix(test_x[, chosen_kmers])

design_matrix_train <- get_interactions_design_3(train_y, train_x_reduced)
design_matrix_test <- get_interactions_design_3(test_y, test_x_reduced)


# LASSO

lasso_cv <- cv.glmnet(design_matrix_train, train_y, family = "binomial")

lambda_min <- lasso_cv$lambda.min

lasso_res <- glmnet(design_matrix_train, train_y, family = "binomial",
                    lambda = lambda_min, thresh = 10^(-10), maxit = 10^(8))

predicted <- predict.glmnet(lasso_res, newx = as.matrix(design_matrix_test), type = "response")
predicted <- as.vector(inv_logit(predicted))

rocobj <- roc(test_y, as.vector(predicted))
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))



# elastic net

elastic_net <- glmnet(design_matrix_train, train_y, family = "binomial", 
                      alpha = 0.5, lambda = lambda_min, thresh = 10^(-10), 
                      maxit = 10^(8))

predicted <- predict.glmnet(elastic_net, newx = as.matrix(design_matrix_test), type = "response")
predicted <- as.vector(apply(predicted, 2, inv_logit))

rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))


# ridge

ridge <- glmnet(design_matrix_train, train_y, family = "binomial", 
                lambda = lambda_min, alpha = 0, thresh = 10^(-10), maxit = 10^(8))

predicted <- predict.glmnet(ridge, newx = as.matrix(design_matrix_test), type = "response")
predicted <- as.vector(apply(predicted, 2, inv_logit))

rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))

## slope

res_slope <- SLOPE(design_matrix_train, train_y, family = "binomial", 
                   lambda = "bh", alpha = 0.01)


predicted <- predict(res_slope, as.matrix(design_matrix_test), type = "response")

rocobj <- roc(test_y, predicted)
ggroc(rocobj)
auc(test_y, predicted)
confusionMatrix(as.factor(as.numeric(predicted > 0.5)), as.factor(test_y))


################# fast forward


dat <- prepare_data(train_y, train_x)

dat_reduced <- dat %>%
  reduce_matrix(minpv = 0.2)

res_maic <- dat_reduced %>%
  fast_forward(crit = "maic")


saveRDS(res_maic[["model"]], "llps_dat/maic_vars.RDS")


res_aic <- dat_reduced %>%
  fast_forward(crit = "aic")

saveRDS(res_aic[["model"]], "llps_dat/aic_vars.RDS")
