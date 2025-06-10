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
Done in r

Filtering outliers:
```bash
plink --bfile gwas_QC --het --out gwas_QC
Filtering done in r
plink --bfile gwas_QC --remove wrong_het.txt --make-bed --out gwas_QC _het 
```

Identifying and removing individuals that are IBD:
```bash
plink --bfile gwas_QC_het --indep-pairwise 500kb 5 0.2 --out gwas_QC_het
plink --bfile gwas_QC_het --extract gwas_QC_het.prune.in --genome --min 0.185 --out gwas_QC_het
filtering done in r
plink --bfile gwas_QC_het --remove wrong_ibd.txt --make-bed --out gwas_QC_ibd
```

Missing height:
Filtering done in r
```bash
plink --bfile  GWA_QC_ibd --remove height_out.txt –allow_no_sex  --make-bed --out gwas_full_height
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
