# GWAS

## software
Plink and R was used for QC, PCA and the LM association model. GCTA used for the LMM association study. R was used for PCA vidualization and LASSO, EM-LASSO and PRS modelling

## QC
Filtering out dubious gender 
```bash
plink --bfile gwas_data --check-sex --out gwas_data
head gwas_data.sexcheck
grep -v "OK" gwas_data.sexcheck > wrong_sex.txt
plink --bfile gwas_data --remove wrong_sex.txt --make-bed --out gwas_QC
```

Filtering out NAs
```r
na_rows <- apply(is.na(d_het),1,any)
rows <- d_het[na_rows,]
d_het <- na.omit(d_het)
```

Filtering outliers:
```bash
plink --bfile gwas_QC --het --out gwas_QC
```
```r
mean_het <- mean(d_het$Het)
sd_het <- sd(d_het$Het)
het_lower <- mean_het - 3*sd_het
het_upper <- mean_het + 3*sd_het
d_m <- d_het %>% filter(Het<het_lower | Het > het_upper)
d_p <- rbind(rows,d_m)
cat("Number of individuals filtered out due to heterozygosity: ", nrow(d_p), "\n")
```
```bash
plink --bfile gwas_QC --remove wrong_het.txt --make-bed --out gwas_QC _het 
```

Identifying and removing individuals that are IBD:
```bash
plink --bfile gwas_QC_het --indep-pairwise 500kb 5 0.2 --out gwas_QC_het
plink --bfile gwas_QC_het --extract gwas_QC_het.prune.in --genome --min 0.185 --out gwas_QC_het
```r
members <- ibd$FID1
members <- unique(members)
cat("Number of individuals filtered out due to IBD: ", length(members), "\n")
```
plink --bfile gwas_QC_het --remove wrong_ibd.txt --make-bed --out gwas_QC_ibd
```

Missing height:
```r
original_n <- nrow(ind)
height_out <- ind %>% anti_join(height, by = c("V1"="V1"))
removed_n <- nrow(height_out)
remaining_n <- original_n - removed_n
cat("Remaining individuals: ", remaining_n, "\n")
```
```bash
plink --bfile  GWA_QC_ibd --remove height_out.txt –allow_no_sex  --make-bed --out gwas_full_height
```

Adding phenotype and genotype data
```r
fam_data<-merge(fam_data,height,by="V1")
fam_data$V6<-fam_data$V2.y
fam_data<-fam_data[,!(names(fam_data)%in% "V2.y")]
colnames(fam_data) <- c("FID", "IID", "PID", "MID", "Sex", "Phenotype")
```

## PCA:
Prune for LD
```bash
plink --bfile gwas_full_height --indep-pairwise 50 5 0.2 --out pca_prune
```

Making the PCA using PLINK:
```bash
plink --bfile GWA_full_height --extract pca_prune.prune.in --pca 10 --out pca_result 
```

Viasualizing PCA
```r
pca <- read.table("pca_result.eigenvec", header = FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:(ncol(pca) - 2)))

fam <- read.table("GWA_full_height.fam")
colnames(fam) <- c("FID", "IID", "PID", "MID", "Sex", "Phenotype")

pca_annot <- merge(pca, fam[, c("FID", "IID", "Sex")], by = c("FID", "IID"))
pca_annot$Sex <- factor(pca_annot$Sex, levels = c(1, 2), labels = c("Male", "Female"))

ggplot(pca_annot, aes(x = PC1, y = PC2, color = Sex)) +
  geom_point(alpha = 0.8, size = 2.5) +
  labs(
    title = "Principal Component Analysis of Genotypic Data",
    subtitle = "Colored by self-reported sex",
    x = "Principal Component 1",
    y = "Principal Component 2",
    color = "Sex"
  ) +
  theme_minimal(base_size = 14)
```

# Making covariates data for GWAS
```r
pca <- read.table("pca_result.eigenvec", header = FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:5))

fam <- read.table("GWA_full_height.fam")
colnames(fam) <- c("FID", "IID", "PID", "MID", "Sex", "Phenotype")

covariates <- merge(pca, fam[, c("FID", "IID", "Sex")], by = c("FID", "IID"))

write.table(covariates, "covariates.txt", quote = FALSE, row.names = FALSE, sep = "\t")

cov <- read.table("covariates.txt", header = TRUE)

# Split
covar_discrete <- cov[, c("FID", "IID", "Sex")]
covar_quant <- cov[, c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5")]

# Save
write.table(covar_discrete, "covar_discrete.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(covar_quant, "covar_quant.txt", sep = "\t", quote = FALSE, row.names = FALSE)
```

## GWAS:
LM GWAS
```bash
plink --bfile GWA_full_height \
   		--covar covariates.txt \
   		--covar-name Sex,PC1,PC2,PC3,PC4,PC5 \
   		--linear \
  		--allow-no-sex \
   		--out gwas_height_corrected
```

Creating ML data:
```bash
plink --bfile GWA_full_height \
      --extract plink_top_snps.txt \
      --recode A \
      --out snp_matrix_ml \
```

# Plotting the LM association

```r
library(data.table)
library(qqman)


assoc <- fread("gwas_height_corrected.assoc.linear")

assoc <- assoc[!is.na(P)]                          
assoc <- assoc[P > 1e-300]                         
assoc[, CHR := as.numeric(CHR)]                    
assoc[, SNP := as.character(SNP)]
assoc <- assoc[CHR %in% 1:22]                      
assoc <- assoc[, .(CHR, BP, SNP, P)]               

assoc[, qval := p.adjust(P, method = "fdr")]
assoc[, FDRsig := qval < 0.05]
assoc[, bonf := P < 0.05 / .N]
assoc[, genomewide := P < 5e-8]

assoc <- assoc[!duplicated(SNP)]
assoc <- assoc[order(P), .SD[1], by = SNP]         # Retain best p-value per SNP

setorder(assoc, CHR, BP)

print(table(cut(assoc$P, breaks = c(0, 1e-100, 1e-50, 1e-10, 0.05, 0.5, 1))))

ylim_cap <- ceiling(max(-log10(assoc$P), na.rm = TRUE)) + 1
col_vec <- rep(c("gray30", "gray60"), length.out = length(unique(assoc$CHR)))

png("manhattan_basic.png", width = 1200, height = 800)
manhattan(assoc,
          chr = "CHR", bp = "BP", snp = "SNP", p = "P",
          genomewideline = -log10(5e-8),
          suggestiveline = -log10(1e-4),
          col = col_vec,
          cex = 0.6,
          ylim = c(0, ylim_cap),
          main = "LM GWAS Manhattan Plot: Height")
dev.off()

png("manhattan_fdr.png", width = 1200, height = 800)
manhattan(assoc,
          chr = "CHR", bp = "BP", snp = "SNP", p = "P",
          highlight = assoc[FDRsig == TRUE, SNP],
          highlight.col = "red3",
          genomewideline = -log10(5e-8),
          suggestiveline = -log10(1e-4),
          col = col_vec,
          cex = 0.6,
          ylim = c(0, ylim_cap),
          main = "GWAS Manhattan Plot with FDR-Adjusted Hits")
dev.off()

png("qqplot.png", width = 800, height = 800)
qq(assoc$P, main = "QQ Plot of GWAS P-values")
dev.off()

cat("Number of SNPs with P < 5e-8:", nrow(assoc[P < 5e-8]), "\n")
cat("Number of SNPs significant after FDR:", nrow(assoc[FDRsig == TRUE]), "\n")
cat("Number of SNPs significant after Bonferroni:", nrow(assoc[bonf == TRUE]), "\n")
```

# Plotting the MLM association
```r
mlm <- fread("gwas_height_mlma.mlma")

mlm <- mlm[!is.na(p) & p > 1e-300]
mlm[, CHR := as.numeric(Chr)]
mlm[, SNP := as.character(SNP)]
mlm[, BP := bp]
mlm <- mlm[CHR %in% 1:22]  # autosomes only
mlm <- mlm[, .(CHR, BP, SNP, P = p)]  # standardize column names

mlm[, qval := p.adjust(P, method = "fdr")]
mlm[, FDRsig := qval < 0.05]

mlm[, bonf := P < 0.05 / .N]
mlm[, genomewide := P < 5e-8]

mlm <- mlm[!duplicated(SNP)]
mlm <- mlm[order(P), .SD[1], by = SNP]

ylim_cap <- 10

png("manhattan_mlm_basic.png", width = 1200, height = 800)
manhattan(mlm,
          chr = "CHR", bp = "BP", snp = "SNP", p = "P",
          genomewideline = -log10(5e-8),
          suggestiveline = -log10(1e-4),
          cex = 0.6,
          ylim = c(0, ylim_cap),
          main = "MLM GWAS Manhattan Plot: Height")
dev.off()

png("manhattan_mlm_fdr.png", width = 1200, height = 800)
manhattan(mlm,
          chr = "CHR", bp = "BP", snp = "SNP", p = "P",
          highlight = mlm[FDRsig == TRUE, SNP],
          genomewideline = -log10(5e-8),
          suggestiveline = -log10(1e-5),
          cex = 1,
          ylim = c(0, ylim_cap),
          main = "MLM Manhattan Plot with FDR Highlighted")
dev.off()

png("qqplot_mlm.png", width = 800, height = 800)
qq(mlm$P, main = "QQ Plot of MLM GWAS P-values")
dev.off()

cat("Number of SNPs with P < 5e-8:", nrow(mlm[P < 5e-8]), "\n")
cat("Number of SNPs significant after FDR:", nrow(mlm[FDRsig == TRUE]), "\n")
cat("Number of SNPs significant after Bonferroni:", nrow(mlm[bonf == TRUE]), "\n")

```

# Full Machine Learning Pipeline with LASSO, Latent LASSO, and Latent RF


```r
##### LM data #####

# Load and clean the data
raw_data_lm <- fread("snp_matrix_ml.raw")
final_data_lm <- raw_data_lm[, !c("FID", "PAT", "MAT"), with = FALSE]
colnames(final_data_lm)[colnames(final_data_lm) == "PHENOTYPE"] <- "y"
final_data_lm$sex <- final_data_lm$SEX - 1

# Loading the PCs
pcs <- fread("pca_result.eigenvec")
colnames(pcs) <- c("FID", "IID", paste0("PC", 1:10))

# making the final data
final_data_lm <- merge(final_data_lm, pcs[, .(IID, PC1, PC2, PC3, PC4, PC5)], by = "IID")
final_data_lm <- final_data_lm %>% dplyr::select(-IID)
  
# Split data
set.seed(1)
split <- initial_split(final_data_lm, prop = 0.8, strata = "y")
train_lm <- training(split)
test_lm <- testing(split)

# Preprocessing
rec <- recipe(y ~ ., data = train_lm) %>% 
  step_impute_knn(all_predictors(), neighbors = 5)
rec_prep <- prep(rec)
train_processed_lm <- bake(rec_prep, new_data = NULL)
test_processed_lm <- bake(rec_prep, new_data = test_lm)

# Prepare matrices
X_train_lm <- as.matrix(train_processed_lm[, colnames(train_processed_lm) != "y"])
y_train_lm <- train_processed_lm$y
X_test_lm <- as.matrix(test_processed_lm[, colnames(test_processed_lm) != "y"])
y_test_lm <- test_processed_lm$y

##### MLM data #####

# Load and clean the data
raw_data_mlm <- fread("mlm_snps_matrix.raw")
final_data_mlm <- raw_data_mlm[, !c("FID", "PAT", "MAT"), with = FALSE]
colnames(final_data_mlm)[colnames(final_data_mlm) == "PHENOTYPE"] <- "y"
final_data_mlm$sex <- final_data_mlm$SEX - 1

# Loading the PCs
pcs <- fread("pca_result.eigenvec")
colnames(pcs) <- c("FID", "IID", paste0("PC", 1:10))

# making the final data
final_data_mlm <- merge(final_data_mlm, pcs[, .(IID, PC1, PC2, PC3, PC4, PC5)], by = "IID")
final_data_mlm <- final_data_mlm %>% dplyr::select(-IID)

# Split data
set.seed(1)
split <- initial_split(final_data_mlm, prop = 0.8, strata = "y")
train_mlm <- training(split)
test_mlm <- testing(split)

# Preprocessing
rec <- recipe(y ~ ., data = train_mlm) %>% 
  step_impute_knn(all_predictors(), neighbors = 5) %>% 
  step_nzv(all_predictors()) %>% 
  step_normalize(all_predictors())
rec_prep <- prep(rec)
train_processed_mlm <- bake(rec_prep, new_data = NULL)
test_processed_mlm <- bake(rec_prep, new_data = test_mlm)

# Prepare matrices
X_train_mlm <- as.matrix(train_processed_mlm[, colnames(train_processed_mlm) != "y"])
y_train_mlm <- train_processed_mlm$y
X_test_mlm <- as.matrix(test_processed_mlm[, colnames(test_processed_mlm) != "y"])
y_test_mlm <- test_processed_mlm$y

# ==============================================================================
# LASSO
# ==============================================================================

##### LM #####

cat("\\nRunning LASSO gor LM GWAS...\\n")
lasso_model_lm <- cv.glmnet(X_train_lm, y_train_lm, alpha = 0.5)
y_pred_lasso_lm <- predict(lasso_model_lm, newx = X_test_lm, s = "lambda.min")
(lasso_r2_lm <- 1 - mean((y_pred_lasso_lm - y_test_lm)^2) / var(y_test_lm))
(lasso_rmse <- sqrt(mean((y_pred_lasso_lm - y_test_lm)^2)) )

##### MLM #####

cat("\\nRunning LASSO gor MLM GWAS...\\n")
lasso_model_mlm <- cv.glmnet(X_train_mlm, y_train_mlm, alpha = 0.5)
y_pred_lasso_mlm <- predict(lasso_model_mlm, newx = X_test_mlm, s = "lambda.min")
(lasso_r2_mlm <- 1 - mean((y_pred_lasso_mlm - y_test_mlm)^2) / var(y_test_mlm))
(lasso_rmse_mlm <- sqrt(mean((y_pred_lasso_mlm - y_test_mlm)^2)) )

# ==============================================================================
# Latent Structure Model
# ==============================================================================

##### LM #####

# Parameters
K <- 3
n <- nrow(X_train_lm)

# Initialization with GMM
gmm_init <- Mclust(X_train_lm[,223:227], G = K)
plot(gmm_init, what = "classification")
z_probs <- gmm_init$z  # This is already soft: n × K matrix with rows summing to 1

# Containers
models <- vector("list", K)
pred_mat <- matrix(NA, nrow = nrow(X_train_lm), ncol = K)
sigma2 <- rep(1, K)

# vectorized loglikelihood function
loglikelihood <- function(y, yhat, sigma2) {
  -0.5 * log(2 * pi * sigma2) - 0.5 * ((y - yhat)^2 / sigma2)
}

# EM loop
for (iter in 1:300) {
  for (k in 1:K) {
    models[[k]] <- cv.glmnet(X_train_lm, y_train_lm, weights = z_probs[, k], alpha = 0)
    pred_mat[, k] <- predict(models[[k]], newx = X_train_lm, s = "lambda.1se")[,1]
    sigma2[k] <- mean((y_train_lm - pred_mat[, k])^2)
  }
  
  log_resp <- matrix(NA, nrow = nrow(X_train_lm), ncol = K)
  for (k in 1:K) {
    log_resp[, k] <- loglikelihood(y_train_lm, pred_mat[, k], sigma2[k])
  }
  
  max_log <- rowMaxs(log_resp)
  log_resp <- log_resp - max_log
  resp <- exp(log_resp)
  z_probs_new <- resp / rowSums(resp)
  
  if (max(abs(z_probs - z_probs_new)) < 1e-6) {
    cat("Converged at iteration", iter, "\\n")
    break
  }
  
  z_probs <- z_probs_new
}

# Cluster Classifier (XGBoost)
cluster_labels_lm <- max.col(z_probs) - 1
dtrain_cluster_lm <- xgb.DMatrix(data = X_train_lm, label = cluster_labels_lm)
cluster_classifier_lm <- xgboost(
  data = dtrain_cluster_lm,
  objective = "multi:softprob",
  num_class = K,
  nrounds = 150,
  max_depth = 4,
  eta = 0.1,
  verbose = 0
)

# Predict soft cluster probabilities on test set
soft_probs_lm <- predict(cluster_classifier_lm, newdata = xgb.DMatrix(X_test_lm))
soft_probs_lm <- matrix(soft_probs_lm, ncol = K, byrow = TRUE)

# Predictions from each expert model
test_preds_lm <- sapply(1:K, function(k) {
  predict(models[[k]], newx = X_test_lm, s = "lambda.min")
})

# Weighted average prediction
latent_pred_soft_lm <- rowSums(test_preds_lm * soft_probs_lm)

# Evaluation
latent_r2_soft_lm <- 1 - mean((latent_pred_soft_lm - y_test_lm)^2) / var(y_test_lm)
latent_rmse_soft_lm <- sqrt(mean((latent_pred_soft_lm - y_test_lm)^2))

cat("Latent Ridge (Soft Assignment) R²:", round(latent_r2_soft_lm, 3), "\n")
cat("Latent Ridge (Soft Assignment) RMSE:", round(latent_rmse_soft_lm, 3), "\n")

##### MLM #####

# Parameters
K <- 3
n <- nrow(X_train_mlm)

# Initialization with GMM
gmm_init <- Mclust(X_train_mlm[,105:109], G = K)
plot(gmm_init, what = "classification")
z_probs <- gmm_init$z  # This is already soft: n × K matrix with rows summing to 1

# Containers
models <- vector("list", K)
pred_mat <- matrix(NA, nrow = nrow(X_train_mlm), ncol = K)
sigma2 <- rep(1, K)

# vectorized loglikelihood function
loglikelihood <- function(y, yhat, sigma2) {
  -0.5 * log(2 * pi * sigma2) - 0.5 * ((y - yhat)^2 / sigma2)
}

# EM loop
for (iter in 1:400) {
  for (k in 1:K) {
    models[[k]] <- cv.glmnet(X_train_mlm, y_train_mlm, weights = z_probs[, k], alpha = 0)
    pred_mat[, k] <- predict(models[[k]], newx = X_train_mlm, s = "lambda.1se")[,1]
    sigma2[k] <- mean((y_train_mlm - pred_mat[, k])^2)
  }
  
  log_resp <- matrix(NA, nrow = nrow(X_train_mlm), ncol = K)
  for (k in 1:K) {
    log_resp[, k] <- loglikelihood(y_train_mlm, pred_mat[, k], sigma2[k])
  }
  
  max_log <- rowMaxs(log_resp)
  log_resp <- log_resp - max_log
  resp <- exp(log_resp)
  z_probs_new <- resp / rowSums(resp)
  
  if (max(abs(z_probs - z_probs_new)) < 1e-6) {
    cat("Converged at iteration", iter, "\\n")
    break
  }
  
  z_probs <- z_probs_new
}

# Cluster Classifier (XGBoost)
cluster_labels_mlm <- max.col(z_probs) - 1
dtrain_cluster_mlm <- xgb.DMatrix(data = X_train_mlm, label = cluster_labels_mlm)
cluster_classifier_mlm <- xgboost(
  data = dtrain_cluster_mlm,
  objective = "multi:softprob",
  num_class = K,
  nrounds = 150,
  max_depth = 4,
  eta = 0.1,
  verbose = 0
)

# Predict soft cluster probabilities on test set
soft_probs_mlm <- predict(cluster_classifier_mlm, newdata = xgb.DMatrix(X_test_mlm))
soft_probs_mlm <- matrix(soft_probs_mlm, ncol = K, byrow = TRUE)

# Predictions from each expert model
test_preds_mlm <- sapply(1:K, function(k) {
  predict(models[[k]], newx = X_test_mlm, s = "lambda.min")
})

# Weighted average prediction
latent_pred_soft_mlm <- rowSums(test_preds_mlm * soft_probs_mlm)

# Evaluation
latent_r2_soft_mlm <- 1 - mean((latent_pred_soft_mlm - y_test_mlm)^2) / var(y_test_mlm)
latent_rmse_soft_mlm <- sqrt(mean((latent_pred_soft_mlm - y_test_mlm)^2))

cat("Latent Ridge (Soft Assignment) R²:", round(latent_r2_soft_mlm, 3), "\n")
cat("Latent Ridge (Soft Assignment) RMSE:", round(latent_rmse_soft_mlm, 3), "\n")

# ==============================================================================
# Making the PRS model score
# ==============================================================================

prs <- fread("prs_pgs001229.profile")
pheno <- fread("height_final.txt")  # Must have FID, IID, height
covars <- fread("covariates.txt")  # Should contain FID, IID, PC1–PC5, Sex

dim(prs)
dim(pheno)
dim(covars)

merged <- merge(prs, pheno, by = c("FID", "IID"))
merged_full <- merge(merged, covars, by = c("FID", "IID"))

set.seed(1)
split <- initial_split(merged_full, prop = 0.8, strata = "PHENO")
train <- training(split)
test <- testing(split)

model_adj <- lm(PHENO ~ SCORE + PC1 + PC2 + PC3 + PC4 + PC5 + Sex, data = train)
summary(model_adj)

y_pred_prs <- predict(model_adj, newdata = test)
(PRS_r2 <- 1 - mean((y_pred_prs - test$PHENO)^2) / var(test$PHENO))
(PRS_rmse <- sqrt(mean((y_pred_prs - test$PHENO)^2)))

# ==============================================================================
# Plotting all performances
# ==============================================================================

# Shared limits across all models
lim_min <- min(c(y_test_lm, y_pred_lasso_lm))
lim_max <- max(c(y_test_lm, y_pred_lasso_lm))

# Set up 2x3 plotting grid
par(mfrow = c(1, 5), mar = c(4, 4, 4, 2))

# LASSO LM
plot(y_test_lm, y_pred_lasso_lm,
     main = "LASSO LM (R² = ", round(lasso_r2_lm, 2), ")",
     xlab = "Actual Height", ylab = "Predicted Height",
     xlim = c(lim_min, lim_max), ylim = c(lim_min, lim_max),
     pch = 19,
     cex = 1.5, 
     cex.main = 2,
     col = rgb(0, 0, 1, 0.3))
abline(0, 1, col = "red", lwd = 2)
grid()

# LASSO MLM
plot(y_test_mlm, y_pred_lasso_mlm,
     main = "LASSO MLM (R² = ", round(lasso_r2_mlm, 2), ")",
     xlab = "Actual Height", ylab = "Predicted Height",
     xlim = c(lim_min, lim_max), ylim = c(lim_min, lim_max),
     pch = 19,
     cex.main = 2,
     cex = 1.5, col = rgb(0, 0.6, 0, 0.3))
abline(0, 1, col = "red", lwd = 2)
grid()

# Latent LM
plot(y_test_lm, latent_pred_soft_lm,
     main = "Latent LM (R² = ", round(latent_r2_soft_lm, 2), ")",
     xlab = "Actual Height", ylab = "Predicted Height",
     xlim = c(lim_min, lim_max), ylim = c(lim_min, lim_max),
     pch = 19,
     cex = 1.5, 
     cex.main = 2,
     col = rgb(0.8, 0.3, 0, 0.3))
abline(0, 1, col = "red", lwd = 2)
grid()

# Latent MLM
plot(y_test_mlm, latent_pred_soft_mlm,
     main = "Latent MLM (R² = ", round(latent_r2_soft_mlm, 2), ")",
     xlab = "Actual Height", ylab = "Predicted Height",
     xlim = c(lim_min, lim_max), ylim = c(lim_min, lim_max),
     pch = 19,
     cex = 1.5, 
     cex.main = 2,
     col = rgb(0.6, 0.2, 0.6, 0.3))
abline(0, 1, col = "red", lwd = 2)
grid()

# PRS Model (external)
plot(test$PHENO, y_pred_prs,
     main = "External PRS (R² = ", round(PRS_r2, 2), ")",
     xlab = "Actual Height", ylab = "Predicted Height",
     xlim = c(lim_min, lim_max), ylim = c(lim_min, lim_max),
     pch = 19,
     cex = 1.5, 
     cex.main = 2,
     col = rgb(0.2, 0.4, 0.8, 0.3))
abline(0, 1, col = "red", lwd = 2)
grid()

# Reset plotting window
par(mfrow = c(1, 1))
```



