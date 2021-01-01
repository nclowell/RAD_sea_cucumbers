## RDA 20201216

#######################################################

library(codep)
library(adespatial)
library(adegraphics)
library(vegan)
library(ape)
library(car)
library(adegenet)
library(tidyverse) 
library(psych)
library(MonteCarlo)
library(data.table)

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

make_zero <- function(adjR2){
  if(is.na(adjR2) == FALSE){
    if(adjR2 < 0){
      return(0)
    } else {
      return(adjR2)
    }
  }
  else {
    return(adjR2)
  }
}

#######################################################
# Read in genetic data, environmental data, and site coordinates

# set wd
setwd("~/Desktop/RDA")

### read in site allele frequencies
afs <- read_csv("PC_pops_noEI_AFs_array.csv", col_names = FALSE) # geographic order
dim(afs) # 8 sites (rows) and 2065 SNPs (columns)
afs[1:8,1:8] # peek
# A tibble: 8 x 8
# X1    X2    X3    X4    X5    X6    X7    X8
# <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#   1 1     0.895 0.895 0.789 0.868 0.941 0.658 1    
# 2 0.923 0.923 0.949 0.833 0.961 0.946 0.718 0.838
# 3 0.95  0.85  0.9   0.8   0.9   1     0.55  0.813
# 4 0.938 0.893 0.9   0.797 0.935 0.96  0.694 0.893
# 5 0.883 0.891 0.936 0.681 0.915 0.944 0.649 0.897
# 6 0.92  0.95  0.917 0.65  0.948 0.885 0.72  0.833
# 7 0.917 0.905 0.929 0.798 0.854 0.91  0.75  0.9  
# 8 0.868 0.947 0.905 0.855 0.878 0.971 0.763 0.821

### hellinger transformation of allele frequencies
afs_hell <- decostand(afs,
                      method="hellinger",
                      MARGIN=1, # SNPs are in rows, so 1 = rows # <--- is this the correct margin?
                      na.rm=TRUE)
afs_hell[1:5,1:5]
# X1         X2         X3         X4         X5
# 1 0.02456706 0.02324153 0.02324153 0.02182185 0.02288827
# 2 0.02363112 0.02363112 0.02396165 0.02244947 0.02411267
# 3 0.02394486 0.02264957 0.02330621 0.02197331 0.02330621
# 4 0.02392601 0.02334504 0.02343636 0.02205454 0.02388772
# 5 0.02312202 0.02322653 0.02380583 0.02030575 0.02353727

### reduce allele frequencies into few PCs using kaiser-guttman criterion (select axes with eigenvalues > mean eigevalue)
geno.pca <- prcomp(afs_hell, scale = T)
screeplot(geno.pca, npcs = 20, type = "lines", bstick = F)
afs_hell_egvals <- geno.pca$sdev^2
KG_afs_hell_egvals <- afs_hell_egvals[afs_hell_egvals > mean(afs_hell_egvals)] 
retain_num <- length(KG_afs_hell_egvals)
geno.pca.axes <- geno.pca$x[,1:retain_num] # 58.29% cumulative genetic variation explained
head(geno.pca.axes)

### read in individual genotype data
gpfile <- read.genepop("parcal_mox001_qc_noEI.gen", ncode = 3L) 
gpfile_dist <- dist(gpfile, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)

### read in environmental data (not yet scaled)
Env <- read_csv("PC_sites_raw_env29.csv")
Env[1:5,1:5]
# A tibble: 5 x 5
# Pop   BO2_salinitymean_ss BO2_ppmean_ss BO2_nitratemean_ss BO2_phosphatemean_ss
# <chr>               <dbl>         <dbl>              <dbl>                <dbl>
#   1 CN_AK               1.11         -1.12               0.926               -0.303
# 2 YT_AK               0.896        -0.968              0.256               -1.15 
# 3 AB_AK               0.698        -0.739              0.297               -1.22 
# 4 BB_BC               0.861        -0.515             -0.689               -1.98 
# 5 SK_WA              -0.491        -0.924             -1.09                 0.123

### reduce env data into few PCs using kaiser-guttman criterion (select axes with eigenvalues > mean eigevalue)
env.pca <- prcomp(Env[,2:ncol(Env)], scale = T)
screeplot(env.pca, npcs = 20, type = "lines", bstick = F)
env_egvals <- env.pca$sdev^2
KG_env_egvals <- env_egvals[env_egvals > mean(env_egvals)] 
retain_num_env <- length(KG_env_egvals)
env.pca.axes <- as.data.frame(env.pca$x[,1:retain_num_env]) # 58.29% cumulative genetic variation explained
head(env.pca.axes)

### which predictors are strongly associated with the first 3 axes?
envPred_load_3PCs <- as_tibble(env.pca$rotation[,1:3], rownames="EnvPred")
#write_csv(envPred_load_3PCs, "EnvPred_loadings_3PCs.csv")
envPred_load_3PCs_abs <- as_tibble(envPred_load_3PCs)
envPred_load_3PCs_abs[2:4] <- abs(envPred_load_3PCs_abs[2:4])

# with the first
envPred_load_PC1_abs <- arrange(envPred_load_3PCs_abs, desc(PC1)) %>%
  select(EnvPred, PC1)

# second
envPred_load_PC2_abs <- arrange(envPred_load_3PCs_abs, desc(PC2)) %>%
  select(EnvPred, PC2)

# third
envPred_load_PC3_abs <- arrange(envPred_load_3PCs_abs, desc(PC3)) %>%
  select(EnvPred, PC3)

### site coords
site_coords <- read_csv("site_coords.csv")
site_coords
# Site    Lat  Long
# <chr> <dbl> <dbl>
#   1 CN_AK  58.2 -152.
# 2 YT_AK  59.7 -140.
# 3 AB_AK  58.4 -135.
# 4 BB_BC  52.2 -128.
# 5 SK_WA  48.3 -124.
# 6 JI_WA  48.5 -123.
# 7 KP_WA  47.7 -123.
# 8 CH_OR  43.3 -124.

# keep just lat and long and save in matrix
Coor <- site_coords %>%
  select(Lat, Long)
Coor_asmat <- as.matrix(Coor)

### get dbMEMs
plot(Coor_asmat[,2], Coor_asmat[,1], asp=1)
DistSpatial<-gcd.hf(Coor_asmat) 
cuke_dbmem <- as.matrix(dbmem(DistSpatial, MEM.autocor = "positive")) # 3

#######################################################
# [RDA1]
# - pop-based
# - partial RDA, condition with dbMEMs
# - all env predictors combined into few PCS
# - all genetic variation in allele frequences

# RESULTS
# VIF: all < 10
# R2 = 0.36 (null = 0.55)
# adjR2 = -0.06 (null = -0004, 0.016 with negatives as zeros)
# quotient real / null adjR2: -3.709452
# difference real - null adjR2: -0.0742034
# ANOVAs insignificant for full model and each axis
# 45 unique SNPs associated with predictors
#######################################################

# run RDA
Y1 <- env.pca.axes
X1 <- afs_hell
Z1 <- cuke_dbmem
rda1 <- rda(X1 ~ . + Condition(Z1), data= Y1)
# Call: rda(formula = X1 ~ PC1 + PC2 + PC3 + Condition(Z1), data = Y1)
# 
# Inertia Proportion Rank
# Total         0.0019313  1.0000000     
# Conditional   0.0004496  0.2328176    1
# Constrained   0.0006924  0.3585424    3
# Unconstrained 0.0007892  0.4086399    3
# Inertia is variance 
# 
# Eigenvalues for constrained axes:
#   RDA1       RDA2       RDA3 
# 0.00027976 0.00024402 0.00016866 
# 
# Eigenvalues for unconstrained axes:
#   PC1       PC2       PC3 
# 0.0003686 0.0002595 0.0001610 

# loadings
coef(rda1)
# RDA1        RDA2        RDA3
# Z1MEM1 -2.6698253 -0.02115795  0.67733077
# Z1MEM2  0.4067097  0.18768327 -0.23355404
# Z1MEM3  0.6215497  0.03201395  0.05607334
# PC1     0.7292418  0.10441969 -0.17913294
# PC2     0.2239738 -0.04209688  0.10870337
# PC3     0.7568632 -0.10539786 -0.23590019

# ANOVA for all axes at once
sig.full.rda1 <- anova.cca(rda1, parallel=getOption("mc.cores")) # default is permutation=999
sig.full.rda1 
# Model: rda(formula = X1 ~ PC1 + PC2 + PC3 + Condition(Z1), data = Y1)
# Df   Variance      F Pr(>F)
# Model     3 0.00069245 0.8774  0.686
# Residual  3 0.00078920

# ANOVA per axis
sig.axis.rda1<- anova.cca(rda1, by="axis", parallel=getOption("mc.cores"))
sig.axis.rda1 
# Model: rda(formula = X1 ~ PC1 + PC2 + PC3 + Condition(Z1), data = Y1)
# Df   Variance      F Pr(>F)
# RDA1      1 0.00027976 1.0635  0.870
# RDA2      1 0.00024402 0.9276  0.849
# RDA3      1 0.00016866 0.6411  0.711
# Residual  3 0.00078920   

vif.cca(rda1)
# Z1       PC1       PC2       PC3 
# 14.005496  7.213699  1.250149  7.541648 

RsquareAdj(rda1)
# $r.squared
# [1] 0.3585424
# 
# $adj.r.squared
# [1] -0.05844713

load.rda1 <- scores(rda1, display="species")
hist(load.rda1[,1], main="Loadings on RDA1")
hist(load.rda1[,2], main="Loadings on RDA2")

# find SNPs
rda1_ax1_cand <- outliers(load.rda1[,1],3)  
length(rda1_ax1_cand) # 22
rda1_ax2_cand  <- outliers(load.rda1[,2],3) 
length(rda1_ax2_cand) # 25
rda1_ax1_cand_df <- cbind.data.frame(rep(1,times=length(rda1_ax1_cand)), names(rda1_ax1_cand), unname(rda1_ax1_cand))
rda1_ax2_cand_df <- cbind.data.frame(rep(2,times=length(rda1_ax2_cand)), names(rda1_ax2_cand), unname(rda1_ax2_cand))
colnames(rda1_ax1_cand_df) <- colnames(rda1_ax2_cand_df) <- c("axis","snp","loading")
rda1_cand <- rbind(rda1_ax1_cand_df, rda1_ax2_cand_df)
rda1_cand$snp <- as.character(rda1_cand$snp)
rda1_unique_snps <- length(unique(rda1_cand$snp)) # 45
rda1_snp_var_pairs <- nrow(rda1_cand) # 47

#######################################################
### Null comparison

num_iter <- 100
null_rda1_results <- tibble(iter=numeric(),R2=numeric(),adjR2=numeric())
for(iteration in seq(num_iter)){
  print(iteration)
  # get random af data, using mean and variance in observed data
  afs_hell_SIM <- as.data.table(afs_hell)[, lapply(.SD, function(x) rnorm(n=8, mean(x), sd(x)))]

  # get random env data, reduce to few PCs using KG criterion
  env29_SIM <- as.data.table(Env[,2:ncol(Env)])[, lapply(.SD, function(x) rnorm(n=8, mean(x), sd(x)))]
  env.sim.pca <- prcomp(env29_SIM, scale = T)
  env.sim_egvals <- env.sim.pca$sdev^2
  KG_env.sim_egvals <- env.sim_egvals[env.sim_egvals > mean(env.sim_egvals)] 
  retain_num_env.sim <- length(KG_env.sim_egvals)
  env.sim.pca.axes <- as.data.frame(env.sim.pca$x[,1:retain_num_env.sim])

  # get random coordinates & dbMEMs
  coors_SIM <- as.data.table(Coor_asmat)[, lapply(.SD, function(x) rnorm(n=8, mean(x), sd(x)))]
  DistSpatial.sim <-gcd.hf(coors_SIM) 
  cuke_dbmem_sim <- as.matrix(dbmem(DistSpatial.sim, MEM.autocor = "positive"))
  
  # run null RDA
  Y1_sim <- env.sim.pca.axes
  X1_sim <- afs_hell_SIM
  Z1_sim <- cuke_dbmem_sim
  rda1_sim <- rda(X1_sim ~ . + Condition(Z1_sim), data= Y1_sim)

  this.R2 <- RsquareAdj(rda1_sim)$r.squared
  this.adjR2 <- RsquareAdj(rda1_sim)$adj.r.squared

  null_rda1_results <- null_rda1_results %>%
    add_row(iter=iteration, R2=this.R2, adjR2=this.adjR2)
}

rda1_sim_mean_R2 <- mean(null_rda1_results$R2) # 0.5457491
rda1_sim_mean_adjR2_wneg <- mean(null_rda1_results$adjR2,na.rm = TRUE) # excluding NAs, keeping neg, -0.004012717
rda1_sim_adjR2_negAs0 <- null_rda1_results %>%
  rowwise() %>%
  mutate(aszero = make_zero(adjR2)) # excluding NAs, keeping neg, -0.004012717
rda1_sim_mean_adjR2_negAs0 <- mean(rda1_sim_adjR2_negAs0$aszero, na.rm = TRUE) # 0.01575627

#######################################################
# [RDA1_nospace]
# - pop-based
# - all env predictors combined into few PCS
# - all genetic variation in allele frequences

# RESULTS
# VIF: all < 10
# R2 = 0.47 (null = 0.54)
# adjR2 = 0.07 (null = -0004, 0.0079 with negatives as zeros)
# quotient real / null adjR2: 8.543331
# difference real - null adjR2: 0.05991744
# ANOVAs insignificant for full model and each axis
# 43 unique SNPs associated with predictors

#######################################################

# run RDA
Y1 <- env.pca.axes
X1 <- afs_hell
rda1_nospace <- rda(X1 ~ ., data= Y1)
# Call: rda(formula = X1 ~ PC1 + PC2 + PC3, data = Y1)
# 
# Inertia Proportion Rank
# Total         0.0019313  1.0000000     
# Constrained   0.0009026  0.4673489    3
# Unconstrained 0.0010287  0.5326511    4
# Inertia is variance 
# 
# Eigenvalues for constrained axes:
#   RDA1      RDA2      RDA3 
# 0.0004793 0.0002543 0.0001690 
# 
# Eigenvalues for unconstrained axes:
#   PC1       PC2       PC3       PC4 
# 0.0004101 0.0002662 0.0001958 0.0001566 

# loadings
coef(rda1_nospace)
# RDA1        RDA2        RDA3
# PC1 0.06351093 -0.09588060  0.04148108
# PC2 0.04951611  0.07646013  0.10091898
# PC3 0.11905858  0.04036146 -0.08899575

# ANOVA for all axes at once
sig.full.rda1_nospace <- anova.cca(rda1_nospace, parallel=getOption("mc.cores")) # default is permutation=999
sig.full.rda1_nospace 
# Model: rda(formula = X1 ~ PC1 + PC2 + PC3, data = Y1)
# Df   Variance      F Pr(>F)
# Model     3 0.00090259 1.1699  0.168
# Residual  4 0.00102870  

# ANOVA per axis
sig.axis.rda1_nospace <- anova.cca(rda1_nospace, by="axis", parallel=getOption("mc.cores"))
sig.axis.rda1_nospace
# Model: rda(formula = X1 ~ PC1 + PC2 + PC3, data = Y1)
# Df   Variance      F Pr(>F)
# RDA1      1 0.00047930 1.8637  0.125
# RDA2      1 0.00025432 0.9889  0.824
# RDA3      1 0.00016897 0.6570  0.765
# Residual  4 0.00102870  

vif.cca(rda1_nospace)
# PC1 PC2 PC3 
# 1   1   1  
RsquareAdj(rda1_nospace)
# $r.squared
# [1] 0.4673489
# 
# $adj.r.squared
# [1] 0.06786054

load.rda1_nospace <- scores(rda1_nospace, display="species")
hist(load.rda1_nospace[,1], main="Loadings on RDA1")
hist(load.rda1_nospace[,2], main="Loadings on RDA2")

# get candidate SNPs
rda1_ns_ax1_cand <- outliers(load.rda1_nospace[,1],3) # 25
length(rda1_ns_ax1_cand)
rda1_ns_ax2_cand  <- outliers(load.rda1_nospace[,2],3) # 20
length(rda1_ns_ax2_cand)
rda1_ns_ax1_cand_df <- cbind.data.frame(rep(1,times=length(rda1_ns_ax1_cand)), names(rda1_ns_ax1_cand), unname(rda1_ns_ax1_cand))
rda1_ns_ax2_cand_df <- cbind.data.frame(rep(2,times=length(rda1_ns_ax2_cand)), names(rda1_ns_ax2_cand), unname(rda1_ns_ax2_cand))
colnames(rda1_ns_ax1_cand_df) <- colnames(rda1_ns_ax2_cand_df) <- c("axis","snp","loading")
rda1_ns_cand <- rbind(rda1_ns_ax1_cand_df, rda1_ns_ax2_cand_df)
rda1_ns_cand$snp <- as.character(rda1_ns_cand$snp)
rda1_ns_unique_snps <- length(unique(rda1_ns_cand$snp)) # 43
rda1_ns_snp_var_pairs <- nrow(rda1_ns_cand) # 45

#######################################################
### Null comparison

num_iter <- 100
null_rda1_ns_results <- tibble(iter=numeric(),R2=numeric(),adjR2=numeric())
for(iteration in seq(num_iter)){
  print(iteration)
  
  # get random af data, using mean and variance in observed data
  afs_hell_SIM <- as.data.table(afs_hell)[, lapply(.SD, function(x) rnorm(n=8, mean(x), sd(x)))]
  
  # get random env data, reduce to few PCs using KG criterion
  env29_SIM <- as.data.table(Env[,2:ncol(Env)])[, lapply(.SD, function(x) rnorm(n=8, mean(x), sd(x)))]
  env.sim.pca <- prcomp(env29_SIM, scale = T)
  env.sim_egvals <- env.sim.pca$sdev^2
  KG_env.sim_egvals <- env.sim_egvals[env.sim_egvals > mean(env.sim_egvals)] 
  retain_num_env.sim <- length(KG_env.sim_egvals)
  env.sim.pca.axes <- as.data.frame(env.sim.pca$x[,1:retain_num_env.sim])

  # run null RDA
  Y1_sim <- env.sim.pca.axes
  X1_sim <- afs_hell_SIM
  rda1_ns_sim <- rda(X1_sim ~ ., data= Y1_sim)
  
  this.R2 <- RsquareAdj(rda1_ns_sim)$r.squared
  this.adjR2 <- RsquareAdj(rda1_ns_sim)$adj.r.squared
  
  null_rda1_ns_results <- null_rda1_ns_results %>%
    add_row(iter=iteration, R2=this.R2, adjR2=this.adjR2)
}

rda1_ns_sim_mean_R2 <- mean(null_rda1_ns_results$R2) # 0.5375313
rda1_ns_sim_mean_adjR2_wneg <- mean(null_rda1_ns_results$adjR2,na.rm = TRUE) # excluding NAs, keeping neg, -0.002623577
rda1_ns_sim_adjR2_negAs0 <- null_rda1_ns_results %>%
  rowwise() %>%
  mutate(aszero = make_zero(adjR2)) 
rda1_ns_sim_mean_adjR2_negAs0 <- mean(rda1_ns_sim_adjR2_negAs0$aszero, na.rm = TRUE) # 0.007943101

#######################################################
# [RDA_ss]
# - pop-based
# - partial RDA, conditioning for space
# - all sea surface env predictors combined into few PCS
# - all genetic variation in allele frequences

# RESULTS
# VIF: 2 > 10 <- too much collinearity
# R2 = 0.25 (null = 0.5)
# adjR2 = 0.16 (null = -0.0011 including negatives, &  0.016 with negatives as zeros)
# quotient real / null adjR2: 9.819719
# difference real - null adjR2: 0.1431604
# ANOVAs insignificant for full model and each axis
# 38 unique SNPs associated with predictors

#######################################################

### sea surface env vars
subset_ss = subset(Env, 
                   select=c(Name,
                            Var.BO2_salinitymean_ss,
                            Var.BO2_ppmean_ss,
                            Var.BO2_nitratemean_ss,
                            Var.BO2_phosphatemean_ss,
                            Var.BO2_dissoxmean_ss,
                            Var.BO_ph,
                            Var.BO_calcite,
                            Var.BO2_curvelmax_ss,
                            Var.BO2_curvelmin_ss,
                            Var.BO2_curvelrange_ss,
                            Var.BO2_curvelmean_ss,
                            Var.BO2_tempmean_ss,
                            Var.BO2_tempmin_ss,
                            Var.BO2_tempmax_ss,
                            Var.BO2_temprange_ss))

# PCA on sea surface environmental predictors, use just few that explain most variation
env_ss_pca <- prcomp(subset_ss[,2:ncol(subset_ss)], scale = TRUE)
env_ss_egvals <- env_ss_pca$sdev^2
KG_env_ss_egvals <- env_ss_egvals[env_ss_egvals > mean(env_ss_egvals)] 
retain_num_env_ss <- length(KG_env_ss_egvals)
env.ss.pca.axes <- as.data.frame(env_ss_pca$x[,1:retain_num_env_ss]) # 58.29% cumulative genetic variation explained
head(env.ss.pca.axes)

### which predictors are strongly associated with the first 2 axes?
envPred_ss_load_2PCs <- as_tibble(env_ss_pca$rotation[,1:2], rownames="EnvPred")
#write_csv(envPred_ss_load_2PCs, "EnvPred_ss_loadings_2PCs.csv")
envPred_ss_load_2PCs_abs <- as_tibble(envPred_ss_load_2PCs)
envPred_ss_load_2PCs_abs[2:3] <- abs(envPred_ss_load_2PCs_abs[2:3])

# with the first
envPred_ss_load_PC1_abs <- arrange(envPred_ss_load_2PCs_abs, desc(PC1)) %>%
  select(EnvPred, PC1)

# second
envPred_ss_load_P2_abs <- arrange(envPred_ss_load_2PCs_abs, desc(PC2)) %>%
  select(EnvPred, PC2)

# run RDA
Y_ss <- env.ss.pca.axes
X_ss <- afs_hell
Z_ss <- cuke_dbmem  
rda6_ss_sp <- rda(X_ss ~ . + Condition(Z_ss), data= Y_ss)
# Call: rda(formula = X_ss ~ PC1 + PC2 + Condition(Z_ss), data = Y_ss)
# 
# Inertia Proportion Rank
# Total         0.0019313  1.0000000     
# Conditional   0.0011350  0.5877095    3
# Constrained   0.0004094  0.2120054    2
# Unconstrained 0.0003868  0.2002852    2
# Inertia is variance 
# 
# Eigenvalues for constrained axes:
#   RDA1       RDA2 
# 0.00021905 0.00019039 
# 
# Eigenvalues for unconstrained axes:
#   PC1        PC2 
# 0.00021967 0.00016714 

coef(rda6_ss_sp)
# RDA1       RDA2
# Z_ssMEM1 -2.2748565 0.10877293
# Z_ssMEM2  0.2682176 0.13490233
# Z_ssMEM3  0.7452107 0.07860537
# PC1       0.5638651 0.13333499
# PC2      -1.0122092 0.12849850

# ANOVA for all axes at once
sig.full.rda6_ss_sp <- anova.cca(rda6_ss_sp, parallel=getOption("mc.cores")) # default is permutation=999
sig.full.rda6_ss_sp 
# Model: rda(formula = X_ss ~ PC1 + PC2 + Condition(Z_ss), data = Y_ss)
# Df   Variance      F Pr(>F)
# Model     2 0.00048608 1.5671  0.161
# Residual  2 0.00031017   

# ANOVA per axis
sig.axis.rda6_ss_sp <- anova.cca(rda6_ss_sp, by="axis", parallel=getOption("mc.cores"))
sig.axis.rda6_ss_sp
# Model: rda(formula = X_ss ~ PC1 + PC2 + Condition(Z_ss), data = Y_ss)
# Df   Variance      F Pr(>F)
# RDA1      1 0.00027598 1.7795  0.286
# RDA2      1 0.00021009 1.3547  0.336
# Residual  2 0.00031017  

vif.cca(rda6_ss_sp)
# Z_ssMEM1  Z_ssMEM2  Z_ssMEM3       PC1       PC2 
# 42.494430  1.721115  5.492142 13.941342 34.766345 

RsquareAdj(rda6_ss_sp)
# $r.squared
# [1] 0.2516859
# 
# $adj.r.squared
# [1] 0.1593922

load.rda6_ss_sp <- scores(rda6_ss_sp, display="species")
hist(load.rda6_ss_sp[,1], main="Loadings on RDA1")
hist(load.rda6_ss_sp[,2], main="Loadings on RDA2")

# find SNPs
rda6_ss_sp_ax1_cand <- outliers(load.rda6_ss_sp[,1],3) # 19
length(rda6_ss_sp_ax1_cand)
rda6_ss_sp_ax2_cand  <- outliers(load.rda6_ss_sp[,2],3) # 19
length(rda6_ss_sp_ax2_cand) 
rda6_ss_sp_ax1_cand_df <- cbind.data.frame(rep(1,times=length(rda6_ss_sp_ax1_cand)), 
                                           names(rda6_ss_sp_ax1_cand), 
                                           unname(rda6_ss_sp_ax1_cand))
rda6_ss_sp_ax2_cand_df <- cbind.data.frame(rep(2,times=length(rda6_ss_sp_ax2_cand)), 
                                           names(rda6_ss_sp_ax2_cand), 
                                           unname(rda6_ss_sp_ax2_cand))
colnames(rda6_ss_sp_ax1_cand_df) <- colnames(rda6_ss_sp_ax2_cand_df) <- c("axis","snp","loading")
rda6_ss_sp_cand <- rbind(rda6_ss_sp_ax1_cand_df, rda6_ss_sp_ax2_cand_df)
rda6_ss_sp_cand$snp <- as.character(rda6_ss_sp_cand$snp)
rda6_ss_sp_unique_snps <- length(unique(rda6_ss_sp_cand$snp)) # 38
rda6_ss_sp_snp_var_pairs <- nrow(rda6_ss_sp_cand) # 38

### plot axes 1 & 2
plot(rda6_ss, type="n", scaling=2)
points(rda6_ss, display="species", pch=20, cex=0.1, col="gray32", scaling=2)           # the SNPs
points(rda6_ss, display="sites", pch=21, cex=2, 
       col=bg, scaling=2, bg=bg) # the wolves
text(rda6_ss, scaling=2, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=levels(site), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
title("Sea surface: NOT conditioning for space")

#######################################################
### Null comparison

num_iter <- 100
null_rda_ss_results <- tibble(iter=numeric(),R2=numeric(),adjR2=numeric())
for(iteration in seq(num_iter)){
  print(iteration)
  
  # get random af data, using mean and variance in observed data
  afs_hell_SIM <- as.data.table(afs_hell)[, lapply(.SD, function(x) rnorm(n=8, mean(x), sd(x)))]
  
  # get random env data, reduce to few PCs using KG criterion
  env_ss_SIM <- as.data.table(subset_ss[,2:ncol(subset_ss)])[, lapply(.SD, function(x) rnorm(n=8, mean(x), sd(x)))]
  env.sim.pca <- prcomp(env_ss_SIM, scale = T)
  env.sim_egvals <- env.sim.pca$sdev^2
  KG_env.sim_egvals <- env.sim_egvals[env.sim_egvals > mean(env.sim_egvals)] 
  retain_num_env.sim <- length(KG_env.sim_egvals)
  env.sim.pca.axes <- as.data.frame(env.sim.pca$x[,1:retain_num_env.sim])
  
  # get random coordinates & dbMEMs
  coors_SIM <- as.data.table(Coor_asmat)[, lapply(.SD, function(x) rnorm(n=8, mean(x), sd(x)))]
  DistSpatial.sim <-gcd.hf(coors_SIM) 
  cuke_dbmem_sim <- as.matrix(dbmem(DistSpatial.sim, MEM.autocor = "positive"))
  
  # run null RDA
  Y_ss_sim <- env.sim.pca.axes
  X_ss_sim <- afs_hell_SIM
  Z_ss_sim <- cuke_dbmem_sim
  rda_ss_sim <- rda(X_ss_sim ~ . + Condition(Z_ss_sim), data= Y_ss_sim)
  
  this.R2 <- RsquareAdj(rda_ss_sim)$r.squared
  this.adjR2 <- RsquareAdj(rda_ss_sim)$adj.r.squared
  
  null_rda_ss_results <- null_rda_ss_results %>%
    add_row(iter=iteration, R2=this.R2, adjR2=this.adjR2)
}

rda_ss_sim_mean_R2 <- mean(null_rda_ss_results$R2) # 0.4964556
rda_ss_sim_mean_adjR2_wneg <- mean(null_rda_ss_results$adjR2,na.rm = TRUE) # excluding NAs, keeping neg, -0.001115473
rda_ss_sim_adjR2_negAs0 <- null_rda_ss_results %>%
  rowwise() %>%
  mutate(aszero = make_zero(adjR2)) 
rda_ss_sim_mean_adjR2_negAs0 <- mean(rda_ss_sim_adjR2_negAs0$aszero, na.rm = TRUE) # 0.01623185

#######################################################
# [RDA_ss_ns]
# - pop-based
# - not conditioning for space
# - all sea surface env predictors combined into few PCS
# - all genetic variation in allele frequences

# RESULTS
# VIF: all < 10
# R2 = 0.36 (null = 0.51)
# adjR2 = 0.11 (null = including negatives,-0.0013 & 0.0062 with negatives as zeros)
# quotient real / null adjR2: 16.85698
# difference real - null adjR2: 0.09897601
# ANOVAs marginally significant for full model and for first axis
# 42 unique SNPs associated with predictors

#######################################################

# run RDA
Y_ss_ns <- env.ss.pca.axes
X_ss_ns <- afs_hell
rda_ss_ns <- rda(X_ss ~ ., data= Y_ss)
# Call: rda(formula = X_ss ~ PC1 + PC2, data = Y_ss)
# 
# Inertia Proportion Rank
# Total         0.0019313  1.0000000     
# Constrained   0.0006969  0.3608699    2
# Unconstrained 0.0012343  0.6391301    5
# Inertia is variance 
# 
# Eigenvalues for constrained axes:
#   RDA1      RDA2 
# 0.0004799 0.0002171 
# 
# Eigenvalues for unconstrained axes:
#   PC1       PC2       PC3       PC4       PC5 
# 0.0004499 0.0002877 0.0001951 0.0001789 0.0001227 

coef(rda_ss_ns)
# RDA1        RDA2
# PC1  0.04127153 -0.14959179
# PC2 -0.16681401 -0.04602304

# ANOVA for all axes at once
sig.full.rda_ss_ns <- anova.cca(rda_ss_ns, parallel=getOption("mc.cores")) # default is permutation=999
sig.full.rda_ss_ns 
# Model: rda(formula = X_ss ~ PC1 + PC2, data = Y_ss)
# Df   Variance      F Pr(>F)  
# Model     2 0.00069694 1.4116  0.063 .
# Residual  5 0.00123434   

# ANOVA per axis
sig.axis.rda_ss_ns <- anova.cca(rda_ss_ns, by="axis", parallel=getOption("mc.cores"))
sig.axis.rda_ss_ns
# Model: rda(formula = X_ss ~ PC1 + PC2, data = Y_ss)
# Df   Variance      F Pr(>F)  
# RDA1      1 0.00047985 1.9438  0.052 .
# RDA2      1 0.00021709 0.8794  0.606  
# Residual  5 0.00123434   

vif.cca(rda_ss_ns)
# PC1 PC2 
# 1   1 

RsquareAdj(rda_ss_ns)
# $r.squared
# [1] 0.3608699
# 
# $adj.r.squared
# [1] 0.1052178

load.rda_ss_ns <- scores(rda_ss_ns, display="species")
hist(load.rda_ss_ns[,1], main="Loadings on RDA1")
hist(load.rda_ss_ns[,2], main="Loadings on RDA2")

# find SNPs
rda_ss_ns_ax1_cand <- outliers(load.rda_ss_ns[,1],3) # 24
length(rda_ss_ns_ax1_cand)
rda_ss_ns_ax2_cand  <- outliers(load.rda_ss_ns[,2],3) # 19
length(rda_ss_ns_ax2_cand) 
rda_ss_ns_ax1_cand_df <- cbind.data.frame(rep(1,times=length(rda_ss_ns_ax1_cand)), 
                                          names(rda_ss_ns_ax1_cand), 
                                          unname(rda_ss_ns_ax1_cand))
rda_ss_ns_ax2_cand_df <- cbind.data.frame(rep(2,times=length(rda_ss_ns_ax2_cand)), 
                                          names(rda_ss_ns_ax2_cand), 
                                          unname(rda_ss_ns_ax2_cand))
colnames(rda_ss_ns_ax1_cand_df) <- colnames(rda_ss_ns_ax2_cand_df) <- c("axis","snp","loading")
rda_ss_ns_cand <- rbind(rda_ss_ns_ax1_cand_df, rda_ss_ns_ax2_cand_df)
rda_ss_ns_cand$snp <- as.character(rda_ss_ns_cand$snp)
rda_ss_ns_unique_snps <- length(unique(rda_ss_ns_cand$snp)) # 42
rda_ss_ns_snp_var_pairs <- nrow(rda_ss_ns_cand) # 43


#######################################################
### Null comparison

num_iter <- 100
null_rda_ss_ns_results <- tibble(iter=numeric(),R2=numeric(),adjR2=numeric())
for(iteration in seq(num_iter)){
  print(iteration)
  
  # get random af data, using mean and variance in observed data
  afs_hell_SIM <- as.data.table(afs_hell)[, lapply(.SD, function(x) rnorm(n=8, mean(x), sd(x)))]
  
  # get random env data, reduce to few PCs using KG criterion
  env_ss_SIM <- as.data.table(subset_ss[,2:ncol(subset_ss)])[, lapply(.SD, function(x) rnorm(n=8, mean(x), sd(x)))]
  env.sim.pca <- prcomp(env_ss_SIM, scale = T)
  env.sim_egvals <- env.sim.pca$sdev^2
  KG_env.sim_egvals <- env.sim_egvals[env.sim_egvals > mean(env.sim_egvals)] 
  retain_num_env.sim <- length(KG_env.sim_egvals)
  env.sim.pca.axes <- as.data.frame(env.sim.pca$x[,1:retain_num_env.sim])
  
  # run null RDA
  Y_ss_sim <- env.sim.pca.axes
  X_ss_sim <- afs_hell_SIM
  rda_ss_ns_sim <- rda(X_ss_sim ~ ., data= Y_ss_sim)
  
  this.R2 <- RsquareAdj(rda_ss_ns_sim)$r.squared
  this.adjR2 <- RsquareAdj(rda_ss_ns_sim)$adj.r.squared
  
  null_rda_ss_ns_results <- null_rda_ss_ns_results %>%
    add_row(iter=iteration, R2=this.R2, adjR2=this.adjR2)
}

rda_ss_ns_sim_mean_R2 <- mean(null_rda_ss_ns_results$R2) # 0.5062532
rda_ss_ns_sim_mean_adjR2_wneg <- mean(null_rda_ss_ns_results$adjR2,na.rm = TRUE) # excluding NAs, keeping neg, -0.001289408
rda_ss_ns_sim_adjR2_negAs0 <- null_rda_ss_ns_results %>%
  rowwise() %>%
  mutate(aszero = make_zero(adjR2)) 
rda_ss_ns_sim_mean_adjR2_negAs0 <- mean(rda_ss_ns_sim_adjR2_negAs0$aszero, na.rm = TRUE) # 0.006241793

#######################################################
# [RDA_bd]
# - pop-based
# - partial RDA, conditioning for space
# - all bottom depth env predictors combined into few PCS
# - all genetic variation in allele frequences

# RESULTS
# VIF: all < 10
# R2 = 0.18 (null = 0.50)
# adjR2 = -0.081 (null = 0.006,  0.02 with negatives as zeros)
# quotient real / null adjR2: -4.124854
# difference real - null adjR2: -0.1010505
# ANOVAs insignificant for full model and each axis
# 34 unique SNPs associated with predictors

#######################################################

### sea surface env vars
subset_bdmean = subset(Env, 
                       select=c(Name,
                                Var.BO2_salinitymean_bdmean,
                                Var.BO2_ppmean_bdmean,
                                Var.BO2_nitratemean_bdmean,
                                Var.BO2_phosphatemean_bdmean,
                                Var.BO2_dissoxmean_bdmean,
                                Var.BO2_curvelmax_bdmean,
                                Var.BO2_curvelmin_bdmean,
                                Var.BO2_curvelrange_bdmean,
                                Var.BO2_curvelmean_bdmean,
                                Var.BO2_tempmean_bdmean,
                                Var.BO2_tempmin_bdmean,
                                Var.BO2_tempmax_bdmean,
                                Var.BO2_temprange_bdmean))

# PCA on sea surface environmental predictors, use just few that explain most variation
env_bd_pca <- prcomp(subset_bdmean[,2:ncol(subset_bdmean)], scale = TRUE)
env_bd_egvals <- env_bd_pca$sdev^2
KG_env_bd_egvals <- env_bd_egvals[env_bd_egvals > mean(env_bd_egvals)] 
retain_num_env_bd <- length(KG_env_bd_egvals)
env.bd.pca.axes <- as.data.frame(env_bd_pca$x[,1:retain_num_env_bd]) # 58.29% cumulative genetic variation explained
head(env.bd.pca.axes)

### which predictors are strongly associated with the first 2 axes?
envPred_bd_load_2PCs <- as_tibble(env_bd_pca$rotation[,1:2], rownames="EnvPred")
#write_csv(envPred_bd_load_2PCs, "EnvPred_bd_loadings_2PCs.csv")
envPred_bd_load_2PCs_abs <- as_tibble(envPred_bd_load_2PCs)
envPred_bd_load_2PCs_abs[2:3] <- abs(envPred_bd_load_2PCs_abs[2:3])

# with the first
envPred_bd_load_PC1_abs <- arrange(envPred_bd_load_2PCs_abs, desc(PC1)) %>%
  select(EnvPred, PC1)

# second
envPred_bd_load_P2_abs <- arrange(envPred_bd_load_2PCs_abs, desc(PC2)) %>%
  select(EnvPred, PC2)

# run RDA
Y_bd <- env.bd.pca.axes
X_bd <- afs_hell
Z_bd <- cuke_dbmem  
rda_bd <- rda(X_bd ~ . + Condition(Z_bd), data= Y_bd)
# Call: rda(formula = X_bd ~ PC1 + PC2 + Condition(Z_bd), data = Y_bd)
# 
# Inertia Proportion Rank
# Total         0.0019313  1.0000000     
# Conditional   0.0011350  0.5877095    3
# Constrained   0.0003532  0.1829073    2
# Unconstrained 0.0004430  0.2293832    2
# Inertia is variance 
# 
# Eigenvalues for constrained axes:
#   RDA1       RDA2 
# 0.00021866 0.00013458 
# 
# Eigenvalues for unconstrained axes:
#   PC1        PC2 
# 0.00026538 0.00017762 

coef(rda_bd)
# RDA1        RDA2
# Z_bdMEM1 -0.53222296  0.08200251
# Z_bdMEM2  0.25082916  0.10887989
# Z_bdMEM3  0.04740792 -0.17843023
# PC1       0.11460416  0.16131727
# PC2      -0.34936822  0.10630594

# ANOVA for all axes at once
sig.full.rda_bd  <- anova.cca(rda_bd, parallel=getOption("mc.cores")) # default is permutation=999
sig.full.rda_bd 
# Model: rda(formula = X_bd ~ PC1 + PC2 + Condition(Z_bd), data = Y_bd)
# Df   Variance      F Pr(>F)
# Model     2 0.00035325 0.7974  0.657
# Residual  2 0.00044300    

# ANOVA per axis
sig.axis.rda_bd <- anova.cca(rda_bd, by="axis", parallel=getOption("mc.cores"))
sig.axis.rda_bd
# Model: rda(formula = X_bd ~ PC1 + PC2 + Condition(Z_bd), data = Y_bd)
# Df   Variance      F Pr(>F)
# RDA1      1 0.00021866 0.9872  0.804
# RDA2      1 0.00013458 0.6076  0.648
# Residual  2 0.00044300  

vif.cca(rda_bd)
# Z_bdMEM1 Z_bdMEM2 Z_bdMEM3      PC1      PC2 
# 3.319885 1.598161 1.272679 1.661606 3.529119 

RsquareAdj(rda_bd)
# $r.squared
# [1] 0.1829073
# 
# $adj.r.squared
# [1] -0.08133275

load.rda_bd <- scores(rda_bd, display="species")
hist(load.rda_bd[,1], main="Loadings on RDA1")
hist(load.rda_bd[,2], main="Loadings on RDA2")

# find SNPs
rda_bd_ax1_cand <- outliers(load.rda_bd[,1],3) # 19
length(rda_bd_ax1_cand)
rda_bd_ax2_cand  <- outliers(load.rda_bd[,2],3) # 19
length(rda_bd_ax2_cand) 
rda_bd_ax1_cand_df <- cbind.data.frame(rep(1,times=length(rda_bd_ax1_cand)), 
                                       names(rda_bd_ax1_cand), 
                                       unname(rda_bd_ax1_cand))
rda_bd_ax2_cand_df <- cbind.data.frame(rep(2,times=length(rda_bd_ax2_cand)), 
                                       names(rda_bd_ax2_cand), 
                                       unname(rda_bd_ax2_cand))
colnames(rda_bd_ax1_cand_df) <- colnames(rda_bd_ax2_cand_df) <- c("axis","snp","loading")
rda_bd_cand <- rbind(rda_bd_ax1_cand_df, rda_bd_ax2_cand_df)
rda_bd_cand$snp <- as.character(rda_bd_cand$snp)
rda_bd_unique_snps <- length(unique(rda_bd_cand$snp)) # 34
rda_bd_snp_var_pairs <- nrow(rda_bd_cand) # 36


#######################################################
### Null comparison

num_iter <- 100
null_rda_bd_results <- tibble(iter=numeric(),R2=numeric(),adjR2=numeric())
for(iteration in seq(num_iter)){
  print(iteration)
  
  # get random af data, using mean and variance in observed data
  afs_hell_SIM <- as.data.table(afs_hell)[, lapply(.SD, function(x) rnorm(n=8, mean(x), sd(x)))]
  
  # get random env data, reduce to few PCs using KG criterion
  env_bd_SIM <- as.data.table(subset_bdmean[,2:ncol(subset_bdmean)])[, lapply(.SD, function(x) rnorm(n=8, mean(x), sd(x)))]
  env.sim.pca <- prcomp(env_bd_SIM, scale = T)
  env.sim_egvals <- env.sim.pca$sdev^2
  KG_env.sim_egvals <- env.sim_egvals[env.sim_egvals > mean(env.sim_egvals)] 
  retain_num_env.sim <- length(KG_env.sim_egvals)
  env.sim.pca.axes <- as.data.frame(env.sim.pca$x[,1:retain_num_env.sim])
  
  # get random coordinates & dbMEMs
  coors_SIM <- as.data.table(Coor_asmat)[, lapply(.SD, function(x) rnorm(n=8, mean(x), sd(x)))]
  DistSpatial.sim <-gcd.hf(coors_SIM) 
  cuke_dbmem_sim <- as.matrix(dbmem(DistSpatial.sim, MEM.autocor = "positive"))
  
  # run null RDA
  Y_bd_sim <- env.sim.pca.axes
  X_bd_sim <- afs_hell_SIM
  Z_bd_sim <- cuke_dbmem_sim
  rda_bd_sim <- rda(X_bd_sim ~ . + Condition(Z_bd_sim), data= Y_bd_sim)
  
  this.R2 <- RsquareAdj(rda_bd_sim)$r.squared
  this.adjR2 <- RsquareAdj(rda_bd_sim)$adj.r.squared
  
  null_rda_bd_results <- null_rda_bd_results %>%
    add_row(iter=iteration, R2=this.R2, adjR2=this.adjR2)
}

rda_bd_sim_mean_R2 <- mean(null_rda_bd_results$R2) # 0.5010866
rda_bd_sim_mean_adjR2_wneg <- mean(null_rda_bd_results$adjR2,na.rm = TRUE) # excluding NAs, keeping neg, 0.006006115
rda_bd_sim_adjR2_negAs0 <- null_rda_bd_results %>%
  rowwise() %>%
  mutate(aszero = make_zero(adjR2)) 
rda_bd_sim_mean_adjR2_negAs0 <- mean(rda_bd_sim_adjR2_negAs0$aszero, na.rm = TRUE) # 0.01971773

#######################################################
# [RDA_bd_ns]
# - pop-based
# - not conditioning for space
# - all bottom depth env predictors combined into few PCS
# - all genetic variation in allele frequences

# RESULTS
# VIF: all < 10
# R2 = 0.29 (null = 0.49)
# adjR2 = -0.0007 (null = 0.003, 0.016 with negatives as zeros)
# quotient real / null adjR2: -0.03374916
# difference real - null adjR2: -0.02038319
# ANOVAs insignificant for full model and each axis
# 39 unique SNPs associated with predictors

#######################################################

# run RDA
Y_bd <- env.bd.pca.axes
X_bd <- afs_hell
rda_bd_ns <- rda(X_bd ~ ., data= Y_bd)
# Call: rda(formula = X_bd ~ PC1 + PC2, data = Y_bd)
# 
# Inertia Proportion Rank
# Total         0.0019313  1.0000000     
# Constrained   0.0005509  0.2852390    2
# Unconstrained 0.0013804  0.7147610    5
# Inertia is variance 
# 
# Eigenvalues for constrained axes:
#   RDA1      RDA2 
# 0.0003463 0.0002046 
# 
# Eigenvalues for unconstrained axes:
#   PC1       PC2       PC3       PC4       PC5 
# 0.0004803 0.0002978 0.0002524 0.0001930 0.0001570 

coef(rda_bd_ns)
# RDA1       RDA2
# PC1  0.02705168 0.15110985
# PC2 -0.19134988 0.03425544

# ANOVA for all axes at once
sig.full.rda_bd_ns  <- anova.cca(rda_bd_ns, parallel=getOption("mc.cores")) # default is permutation=999
sig.full.rda_bd_ns 
# Model: rda(formula = X_bd ~ PC1 + PC2, data = Y_bd)
# Df   Variance      F Pr(>F)
# Model     2 0.00055088 0.9977  0.465
# Residual  5 0.00138041    

# ANOVA per axis
sig.axis.rda_bd_ns <- anova.cca(rda_bd_ns, by="axis", parallel=getOption("mc.cores"))
sig.axis.rda_bd_ns
# Model: rda(formula = X_bd ~ PC1 + PC2, data = Y_bd)
# Df   Variance      F Pr(>F)
# RDA1      1 0.00034632 1.2544  0.435
# RDA2      1 0.00020456 0.7409  0.821
# Residual  5 0.00138041  

vif.cca(rda_bd_ns)
# PC1 PC2 
# 1   1 

RsquareAdj(rda_bd_ns)
# $r.squared
# [1] 0.285239
# 
# $adj.r.squared
# [1] -0.0006654568

load.rda_bd_ns <- scores(rda_bd_ns, display="species")
hist(load.rda_bd_ns[,1], main="Loadings on RDA1")
hist(load.rda_bd_ns[,2], main="Loadings on RDA2")

# find SNPs
rda_bd_ns_ax1_cand <- outliers(load.rda_bd_ns[,1],3) # 25
length(rda_bd_ns_ax1_cand)
rda_bd_ns_ax2_cand  <- outliers(load.rda_bd_ns[,2],3) # 16
length(rda_bd_ns_ax2_cand) 
rda_bd_ns_ax1_cand_df <- cbind.data.frame(rep(1,times=length(rda_bd_ns_ax1_cand)), 
                                           names(rda_bd_ns_ax1_cand), 
                                           unname(rda_bd_ns_ax1_cand))
rda_bd_ns_ax2_cand_df <- cbind.data.frame(rep(2,times=length(rda_bd_ns_ax2_cand)), 
                                           names(rda_bd_ns_ax2_cand), 
                                           unname(rda_bd_ns_ax2_cand))
colnames(rda_bd_ns_ax1_cand_df) <- colnames(rda_bd_ns_ax2_cand_df) <- c("axis","snp","loading")
rda_bd_ns_cand <- rbind(rda_bd_ns_ax1_cand_df, rda_bd_ns_ax2_cand_df)
rda_bd_ns_cand$snp <- as.character(rda_bd_ns_cand$snp)
rda_bd_ns_unique_snps <- length(unique(rda_bd_ns_cand$snp)) # 39
rda_bd_ns_snp_var_pairs <- nrow(rda_bd_ns_cand) # 41

#######################################################
### Null comparison

num_iter <- 100
null_rda_bd_ns_results <- tibble(iter=numeric(),R2=numeric(),adjR2=numeric())
for(iteration in seq(num_iter)){
  print(iteration)
  
  # get random af data, using mean and variance in observed data
  afs_hell_SIM <- as.data.table(afs_hell)[, lapply(.SD, function(x) rnorm(n=8, mean(x), sd(x)))]
  
  # get random env data, reduce to few PCs using KG criterion
  env_bd_SIM <- as.data.table(subset_bdmean[,2:ncol(subset_bdmean)])[, lapply(.SD, function(x) rnorm(n=8, mean(x), sd(x)))]
  env.sim.pca <- prcomp(env_bd_SIM, scale = T)
  env.sim_egvals <- env.sim.pca$sdev^2
  KG_env.sim_egvals <- env.sim_egvals[env.sim_egvals > mean(env.sim_egvals)] 
  retain_num_env.sim <- length(KG_env.sim_egvals)
  env.sim.pca.axes <- as.data.frame(env.sim.pca$x[,1:retain_num_env.sim])
  
  # run null RDA
  Y_bd_sim <- env.sim.pca.axes
  X_bd_sim <- afs_hell_SIM
  rda_bd_ns_sim <- rda(X_bd_sim ~ . + Condition(Z_bd_sim), data= Y_bd_sim)
  
  this.R2 <- RsquareAdj(rda_bd_ns_sim)$r.squared
  this.adjR2 <- RsquareAdj(rda_bd_ns_sim)$adj.r.squared
  
  null_rda_bd_ns_results <- null_rda_bd_ns_results %>%
    add_row(iter=iteration, R2=this.R2, adjR2=this.adjR2)
}

rda_bd_ns_sim_mean_R2 <- mean(null_rda_bd_ns_results$R2) # 0.4925438
rda_bd_ns_sim_mean_adjR2_wneg <- mean(null_rda_bd_ns_results$adjR2,na.rm = TRUE) # excluding NAs, keeping neg, 0.003053486
rda_bd_ns_sim_adjR2_negAs0 <- null_rda_bd_ns_results %>%
  rowwise() %>%
  mutate(aszero = make_zero(adjR2)) 
rda_bd_ns_sim_mean_adjR2_negAs0 <- mean(rda_bd_ns_sim_adjR2_negAs0$aszero, na.rm = TRUE) # 0.01604612

#######################################################
# [RDA_cv_temp] 
# - pop-based
# - partial RDA, conditioning for space
# - cv and temp predictors into few PCs
# - all genetic variation in allele frequences

# RESULTS
# VIF: 2 > 10 -> too much collinearity
# R2 =  0.32 (null = 0.52)
# adjR2 = 0.087 (null = -0.00013, 0.016 with negatives as zeros)
# quotient real / null adjR2: 5.288327
# difference real - null adjR2: 0.07039983
# ANOVAs insignificant for full model and each axis
# 40 unique SNPs associated with predictors

#######################################################

# note: i tried combining curvel and temp env vars separately in different PCAs, but they are
# SO correlated that the vif were huge

### curvel and temp env vars
subset_cv_temp = subset(Env, 
                       select=c(Name,
                                Var.BO2_curvelmax_bdmean,
                                Var.BO2_curvelmin_bdmean,
                                Var.BO2_curvelrange_bdmean,
                                Var.BO2_curvelmean_bdmean,
                                Var.BO2_curvelmax_ss,
                                Var.BO2_curvelmin_ss,
                                Var.BO2_curvelrange_ss,
                                Var.BO2_curvelmean_ss,
                                Var.BO2_tempmean_bdmean,
                                Var.BO2_tempmin_bdmean,
                                Var.BO2_tempmax_bdmean,
                                Var.BO2_temprange_bdmean,
                                Var.BO2_tempmean_ss,
                                Var.BO2_tempmin_ss,
                                Var.BO2_tempmax_ss,
                                Var.BO2_temprange_ss))


# PCA on cv_temp environmental predictors, use just few that explain most variation
env_cv_temp_pca <- prcomp(subset_cv_temp[,2:ncol(subset_cv_temp)], scale = TRUE)
env_cv_temp_egvals <- env_cv_temp_pca$sdev^2
KG_env_cv_temp_egvals <- env_cv_temp_egvals[env_cv_temp_egvals > mean(env_cv_temp_egvals)] 
retain_num_env_cv_temp <- length(KG_env_cv_temp_egvals)
env.cv_temp.pca.axes <- as.data.frame(env_cv_temp_pca$x[,1:retain_num_env_cv_temp]) # 58.29% cumulative genetic variation explained
head(env.cv_temp.pca.axes)

### which predictors are strongly associated with the first 3 axes?
envPred_cv_temp_load_3PCs <- as_tibble(env_cv_temp_pca$rotation[,1:3], rownames="EnvPred")
write_csv(envPred_cv_temp_load_3PCs, "EnvPred_cv_temp_loadings_3PCs.csv")
envPred_cv_temp_load_3PCs_abs <- as_tibble(envPred_cv_temp_load_3PCs)
envPred_cv_temp_load_3PCs_abs[2:4] <- abs(envPred_cv_temp_load_3PCs_abs[2:4])

# with the first
envPred_cv_temp_load_PC1_abs <- arrange(envPred_cv_temp_load_3PCs_abs, desc(PC1)) %>%
  select(EnvPred, PC1)

# second
envPred_cv_temp_load_P2_abs <- arrange(envPred_cv_temp_load_3PCs_abs, desc(PC2)) %>%
  select(EnvPred, PC2)

# run RDA
Y_cv_temp <- env.cv_temp.pca.axes
X_cv_temp <- afs_hell
Z_cv_temp <- cuke_dbmem  
rda_cv_temp <- rda(X_cv_temp ~ . + Condition(Z_cv_temp), data= Y_cv_temp)
# Call: rda(formula = X_cv_temp ~ PC1 + PC2 + PC3 + Condition(Z_cv_temp), data = Y_cv_temp)
# 
# Inertia Proportion Rank
# Total         0.0019313  1.0000000     
# Conditional   0.0011350  0.5877095    3
# Constrained   0.0006211  0.3216203    3
# Unconstrained 0.0001751  0.0906703    1
# Inertia is variance 
# 
# Eigenvalues for constrained axes:
#   RDA1       RDA2       RDA3 
# 0.00027659 0.00021678 0.00012777 
# 
# Eigenvalues for unconstrained axes:
#   PC1 
# 0.00017511 

coef(rda_cv_temp)
# RDA1          RDA2       RDA3
# Z_cv_tempMEM1 -1.80167123 -0.1031260291  0.9873795
# Z_cv_tempMEM2 -0.08209399  0.2642086070 -0.1359105
# Z_cv_tempMEM3  0.47742066 -0.0005681642  0.1889201
# PC1            0.48905797  0.1395716128 -0.2309374
# PC2            0.73339051 -0.1564404042 -0.1430879
# PC3           -0.49890222  0.0303252446  0.6036349

# ANOVA for all axes at once
sig.full.rda_cv_temp  <- anova.cca(rda_cv_temp, parallel=getOption("mc.cores")) # default is permutation=999
sig.full.rda_cv_temp
# Model: rda(formula = X_cv_temp ~ PC1 + PC2 + PC3 + Condition(Z_cv_temp), data = Y_cv_temp)
# Df   Variance      F Pr(>F)
# Model     3 0.00062114 1.1824  0.367
# Residual  1 0.00017511  

# ANOVA per axis
sig.axis.rda_curvel <- anova.cca(rda_curvel, by="axis", parallel=getOption("mc.cores"))
sig.axis.rda_curvel
# Model: rda(formula = X_curvel ~ PC1 + PC2 + Condition(Z_curvel), data = Y_curvel)
# Df   Variance      F Pr(>F)
# RDA1      1 0.00024979 1.3347  0.572
# RDA2      1 0.00017214 0.9198  0.529
# Residual  2 0.00037432 

vif.cca(rda_cv_temp)
# Z_cv_tempMEM1 Z_cv_tempMEM2 Z_cv_tempMEM3           PC1           PC2           PC3 
# 34.852579      1.760138      3.108973     17.958920     12.090733      9.672038 

RsquareAdj(rda_cv_temp)
# $r.squared
# [1] 0.3216203
# 
# $adj.r.squared
# [1] 0.08681645

load.rda_cv_temp <- scores(rda_cv_temp, display="species")
hist(load.rda_cv_temp[,1], main="Loadings on RDA1")
hist(load.rda_cv_temp[,2], main="Loadings on RDA2")

# find SNPs
rda_cv_temp_ax1_cand <- outliers(load.rda_cv_temp[,1],3) # 21
length(rda_cv_temp_ax1_cand)
rda_cv_temp_ax2_cand  <- outliers(load.rda_cv_temp[,2],3) # 19
length(rda_cv_temp_ax2_cand) 
rda_cv_temp_ax1_cand_df <- cbind.data.frame(rep(1,times=length(rda_cv_temp_ax1_cand)), 
                                            names(rda_cv_temp_ax1_cand), 
                                            unname(rda_cv_temp_ax1_cand))
rda_cv_temp_ax2_cand_df <- cbind.data.frame(rep(2,times=length(rda_cv_temp_ax2_cand)), 
                                            names(rda_cv_temp_ax2_cand), 
                                            unname(rda_cv_temp_ax2_cand))
colnames(rda_cv_temp_ax1_cand_df) <- colnames(rda_cv_temp_ax2_cand_df) <- c("axis","snp","loading")
rda_cv_temp_cand <- rbind(rda_cv_temp_ax1_cand_df, rda_cv_temp_ax2_cand_df)
rda_cv_temp_unique_snps <- length(unique(rda_cv_temp_cand$snp)) # 40
rda_cv_temp_snp_var_pairs <- nrow(rda_cv_temp_cand) # 40

#######################################################
### Null comparison

num_iter <- 100
null_rda_cv_temp_results <- tibble(iter=numeric(),R2=numeric(),adjR2=numeric())
for(iteration in seq(num_iter)){
  print(iteration)
  
  # get random af data, using mean and variance in observed data
  afs_hell_SIM <- as.data.table(afs_hell)[, lapply(.SD, function(x) rnorm(n=8, mean(x), sd(x)))]
  
  # get random env data, reduce to few PCs using KG criterion
  env_cv_temp_SIM <- as.data.table(subset_cv_temp[,2:ncol(subset_cv_temp)])[, lapply(.SD, function(x) rnorm(n=8, mean(x), sd(x)))]
  env.sim.pca <- prcomp(env_cv_temp_SIM, scale = T)
  env.sim_egvals <- env.sim.pca$sdev^2
  KG_env.sim_egvals <- env.sim_egvals[env.sim_egvals > mean(env.sim_egvals)] 
  retain_num_env.sim <- length(KG_env.sim_egvals)
  env.sim.pca.axes <- as.data.frame(env.sim.pca$x[,1:retain_num_env.sim])
  
  # get random coordinates & dbMEMs
  coors_SIM <- as.data.table(Coor_asmat)[, lapply(.SD, function(x) rnorm(n=8, mean(x), sd(x)))]
  DistSpatial.sim <-gcd.hf(coors_SIM) 
  cuke_dbmem_sim <- as.matrix(dbmem(DistSpatial.sim, MEM.autocor = "positive"))
  
  # run null RDA
  Y_cv_temp_sim <- env.sim.pca.axes
  X_cv_temp_sim <- afs_hell_SIM
  Z_cv_temp_sim <- cuke_dbmem_sim
  rda_cv_temp_sim <- rda(X_cv_temp_sim ~ . + Condition(Z_cv_temp_sim), data= Y_cv_temp_sim)
  
  this.R2 <- RsquareAdj(rda_cv_temp_sim)$r.squared
  this.adjR2 <- RsquareAdj(rda_cv_temp_sim)$adj.r.squared
  
  null_rda_cv_temp_results <- null_rda_cv_temp_results %>%
    add_row(iter=iteration, R2=this.R2, adjR2=this.adjR2)
}

rda_cv_temp_sim_mean_R2 <- mean(null_rda_cv_temp_results$R2) # 0.515088
rda_cv_temp_sim_mean_adjR2_wneg <- mean(null_rda_cv_temp_results$adjR2,na.rm = TRUE) # excluding NAs, keeping neg, -0.0001313236
rda_cv_temp_sim_adjR2_negAs0 <- null_rda_cv_temp_results %>%
  rowwise() %>%
  mutate(aszero = make_zero(adjR2)) 
rda_cv_temp_sim_mean_adjR2_negAs0 <- mean(rda_cv_temp_sim_adjR2_negAs0$aszero, na.rm = TRUE) # 0.01641662

#######################################################
# [RDA_cv_temp_ns] 
# - pop-based
# - NOT conditioning for space
# - cv and temp predictors into few PCs
# - all genetic variation in allele frequences

# RESULTS
# VIF: all < 10
# R2 = 0.51 (null = 0.51)
# adjR2 = 0.14 (null = 0.0011, 0.018 with negatives as zeros)
# quotient real / null adjR2: 8.043955
# difference real - null adjR2: 0.1233377
# ANOVAs significant for full model and first axis
# 47 unique SNPs associated with predictors

#######################################################

# run RDA
Y_cv_temp <- env.cv_temp.pca.axes
X_cv_temp <- afs_hell
rda_cv_temp_ns <- rda(X_cv_temp ~ ., data= Y_cv_temp)
# Call: rda(formula = X_cv_temp ~ PC1 + PC2 + PC3, data = Y_cv_temp)
# 
# Inertia Proportion Rank
# Total         0.0019313  1.0000000     
# Constrained   0.0009831  0.5090556    3
# Unconstrained 0.0009482  0.4909444    4
# Inertia is variance 
# 
# Eigenvalues for constrained axes:
#   RDA1      RDA2      RDA3 
# 0.0005186 0.0002780 0.0001865 
# 
# Eigenvalues for unconstrained axes:
#   PC1       PC2       PC3       PC4 
# 0.0003667 0.0002513 0.0001867 0.0001434 

coef(rda_cv_temp_ns)
# RDA1        RDA2       RDA3
# PC1  0.06005018 -0.09777649 0.06485493
# PC2  0.17792936  0.12609886 0.02536143
# PC3 -0.09281092  0.08722682 0.21743982

# ANOVA for all axes at once
sig.full.rda_cv_temp_ns  <- anova.cca(rda_cv_temp_ns, parallel=getOption("mc.cores")) # default is permutation=999
sig.full.rda_cv_temp_ns
# Model: rda(formula = X_cv_temp ~ PC1 + PC2 + PC3, data = Y_cv_temp)
# Df   Variance      F Pr(>F)  
# Model     3 0.00098313 1.3825  0.034 *
#   Residual  4 0.00094815          

# ANOVA per axis
sig.axis.rda_cv_temp_ns <- anova.cca(rda_cv_temp_ns, by="axis", parallel=getOption("mc.cores"))
sig.axis.rda_cv_temp_ns
# Model: rda(formula = X_cv_temp ~ PC1 + PC2 + PC3, data = Y_cv_temp)
# Df   Variance      F Pr(>F)  
# RDA1      1 0.00051858 2.1877  0.035 *
#   RDA2      1 0.00027804 1.1730  0.624  
# RDA3      1 0.00018651 0.7868  0.695  
# Residual  4 0.00094815 

vif.cca(rda_cv_temp_ns)
# PC1 PC2 PC3 
# 1   1   1 

RsquareAdj(rda_cv_temp_ns)
# $r.squared
# [1] 0.5090556
# 
# $adj.r.squared
# [1] 0.1408474

load.rda_cv_temp_ns <- scores(rda_cv_temp_ns, display="species")
hist(load.rda_cv_temp_ns[,1], main="Loadings on RDA1")
hist(load.rda_cv_temp_ns[,2], main="Loadings on RDA2")

# find SNPs
rda_cv_temp_ns_ax1_cand <- outliers(load.rda_cv_temp_ns[,1],3) # 29
length(rda_cv_temp_ns_ax1_cand)
rda_cv_temp_ns_ax2_cand  <- outliers(load.rda_cv_temp_ns[,2],3) # 21
length(rda_cv_temp_ns_ax2_cand) 
rda_cv_temp_ns_ax1_cand_df <- cbind.data.frame(rep(1,times=length(rda_cv_temp_ns_ax1_cand)), 
                                               names(rda_cv_temp_ns_ax1_cand), 
                                               unname(rda_cv_temp_ns_ax1_cand))
rda_cv_temp_ns_ax2_cand_df <- cbind.data.frame(rep(2,times=length(rda_cv_temp_ns_ax2_cand)), 
                                               names(rda_cv_temp_ns_ax2_cand), 
                                               unname(rda_cv_temp_ns_ax2_cand))
colnames(rda_cv_temp_ns_ax1_cand_df) <- colnames(rda_cv_temp_ns_ax2_cand_df) <- c("axis","snp","loading")
rda_cv_temp_ns_cand <- rbind(rda_cv_temp_ns_ax1_cand_df, rda_cv_temp_ns_ax2_cand_df)
rda_cv_temp_ns_unique_snps <- length(unique(rda_cv_temp_ns_cand$snp)) # 47
rda_cv_temp_ns_snp_var_pairs <- nrow(rda_cv_temp_ns_cand) # 50

#######################################################
### Null comparison

num_iter <- 100
null_rda_cv_temp_ns_results <- tibble(iter=numeric(),R2=numeric(),adjR2=numeric())
for(iteration in seq(num_iter)){
  print(iteration)
  
  # get random af data, using mean and variance in observed data
  afs_hell_SIM <- as.data.table(afs_hell)[, lapply(.SD, function(x) rnorm(n=8, mean(x), sd(x)))]
  
  # get random env data, reduce to few PCs using KG criterion
  env_cv_temp_ns_SIM <- as.data.table(subset_cv_temp[,2:ncol(subset_cv_temp)])[, lapply(.SD, function(x) rnorm(n=8, mean(x), sd(x)))]
  env.sim.pca <- prcomp(env_cv_temp_ns_SIM, scale = T)
  env.sim_egvals <- env.sim.pca$sdev^2
  KG_env.sim_egvals <- env.sim_egvals[env.sim_egvals > mean(env.sim_egvals)] 
  retain_num_env.sim <- length(KG_env.sim_egvals)
  env.sim.pca.axes <- as.data.frame(env.sim.pca$x[,1:retain_num_env.sim])
  
  # get random coordinates & dbMEMs
  coors_SIM <- as.data.table(Coor_asmat)[, lapply(.SD, function(x) rnorm(n=8, mean(x), sd(x)))]
  DistSpatial.sim <-gcd.hf(coors_SIM) 
  cuke_dbmem_sim <- as.matrix(dbmem(DistSpatial.sim, MEM.autocor = "positive"))
  
  # run null RDA
  Y_cv_temp_ns_sim <- env.sim.pca.axes
  X_cv_temp_ns_sim <- afs_hell_SIM
  Z_cv_temp_ns_sim <- cuke_dbmem_sim
  rda_cv_temp_ns_sim <- rda(X_cv_temp_ns_sim ~ . + Condition(Z_cv_temp_ns_sim), data= Y_cv_temp_ns_sim)
  
  this.R2 <- RsquareAdj(rda_cv_temp_ns_sim)$r.squared
  this.adjR2 <- RsquareAdj(rda_cv_temp_ns_sim)$adj.r.squared
  
  null_rda_cv_temp_ns_results <- null_rda_cv_temp_ns_results %>%
    add_row(iter=iteration, R2=this.R2, adjR2=this.adjR2)
}

rda_cv_temp_ns_sim_mean_R2 <- mean(null_rda_cv_temp_ns_results$R2) # 0.5066517
rda_cv_temp_ns_sim_mean_adjR2_wneg <- mean(null_rda_cv_temp_ns_results$adjR2,na.rm = TRUE) # excluding NAs, keeping neg, 0.001061491
rda_cv_temp_ns_sim_adjR2_negAs0 <- null_rda_cv_temp_ns_results %>%
  rowwise() %>%
  mutate(aszero = make_zero(adjR2)) 
rda_cv_temp_ns_sim_mean_adjR2_negAs0 <- mean(rda_cv_temp_ns_sim_adjR2_negAs0$aszero, na.rm = TRUE) # 0.01750972



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### comparing SNPs across all RDAs here

rda1_snps <- as.vector(unique(rda1_cand$snp))
rda1_ns_snps <-as.vector(unique(rda1_ns_cand$snp))

rda_ss_snps <- as.vector(unique(rda6_ss_sp_cand$snp))
rda_ss_ns_snps <- as.vector(unique(rda_ss_ns_cand$snp))

rda_bd_snps <- as.vector(unique(rda_bd_cand$snp))
rda_bd_ns_snps <- as.vector(unique(rda_bd_ns_cand$snp))

rda_cv_temp_snps <- as.vector(unique(rda_cv_temp_cand$snp))
rda_cv_temp_ns_snps <- as.vector(unique(rda_cv_temp_ns_cand$snp))

library(VennDiagram)

# all partial RDAs
venn.diagram(
  x = list(rda1_snps, rda_ss_snps, rda_bd_snps, rda_cv_temp_snps ),
  category.names = c("all" , "ss", "bd", "cv_temp"),
  filename = 'vennD_snps_pRDAs.png',
  output=TRUE
)

# all reg RDAS
venn.diagram(
  x = list(rda1_ns_snps, rda_ss_ns_snps, rda_bd_ns_snps, rda_cv_temp_ns_snps ),
  category.names = c("all" , "ss", "bd", "cv_temp"),
  filename = 'vennD_snps_regRDAs.png',
  output=TRUE
)

### % sim between each partial and reg

# all29
length(intersect(rda1_snps,rda1_ns_snps)) # 23 shared

# ss
length(intersect(rda_ss_snps,rda_ss_ns_snps)) # 10 shared

# bd
length(intersect(rda_bd_snps,rda_bd_ns_snps)) # 8 shared

# cv temp
length(intersect(rda_cv_temp_snps,rda_cv_temp_ns_snps)) # 13


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### output file for gene annotation

store_snp_info <- function(rda_name, rda_cand, gene_annotation_tibble){
  rows <- nrow(rda_cand)
  for(j in seq(rows)){
    this_axis <- rda_cand[j,1]
    snp_index_x <- rda_cand[j,2]
    snp_index <- as.numeric(substr(snp_index_x, 2, nchar(snp_index_x)))
    this_snp <- snp_names_2065$V1[snp_index]
    this_snp_loading <- rda_cand[j,3]
    for_gene_annotation <- for_gene_annotation %>%
      add_row(RDA=rda_name, axis=this_axis, snp=this_snp, snp_loading=this_snp_loading)
  }
  return(for_gene_annotation)
}

# load SNP names in order
snp_names_2065 <- read.delim2("PC_2065_biallelic_SNPs.txt", header = FALSE)

# init tibble & store with candidate snps
for_gene_annotation <- tibble(RDA=character(), axis=numeric(), snp=character(), snp_loading=numeric())
rda_names <- c("rda1", "rda1_ns", "rda_ss", "rda_ss_ns", "rda_bd", "rda_bd_ns", "rda_cv_temp", "rda_cv_temp_ns")
rda_cands <- list(rda1_cand, rda1_ns_cand, rda6_ss_sp_cand, rda_ss_ns_cand, rda_bd_cand, rda_bd_ns_cand, rda_cv_temp_cand, rda_cv_temp_ns_cand)
num_rdas <- length(rda_names)
for(i in seq(num_rdas)){
  rda_name <- rda_names[i]
  rda_cand <- rda_cands[[i]]
  for_gene_annotation <- store_snp_info(rda_name, rda_cand, for_gene_annotation)
}
for_gene_annotation 
write_csv(for_gene_annotation, "RDA_SNPs_for_annotation.csv")



