###################### USING OUTFLANK TO IDENTIFY OUTLIERS ##########################
# 20191024 NL
# PURPOSE: Runs OutFLANK 
# INPUTS, made with custom python script convert_genepop_for_outflank.py:
# - SNP matrix file
# - populations file
# - locus names file
#####################################################################################

# set working directory to folder with OutFLANK input files
setwd("E:/dDocent_for_mox/parcal_wd/parcal_mox001/outliers/")

# # install necessary packages
# install.packages("devtools")
# source("http://bioconductor.org/biocLite.R")
# biocLite("qvalue")
# install_github("whitlock/OutFLANK", force = TRUE)
# source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/Fst%20Diploids.R")

# load packages
library(OutFLANK)
library(tidyverse)

# load input files
loci <- read.table("parcal_mox001_loci.txt", header=F)
pops <- read.table("parcal_mox001_pops.txt", header=F)
data <- read.table("parcal_mox001_snpmat.txt", header=F)

corename <- "primSNPs_noINDL_parcal_mox001_md70_maf05_minQ20_mDP10_inames_noreps_HWE_oneSNPhiMAF_gorder"

# read sNP table in as matrix
datamat = as.matrix(data)

# get number of populations
num_pops <- length(unique(pops$V1))

# make Fst DataFrame
FstDataFrame <- MakeDiploidFSTMat(SNPmat = datamat, locusNames = loci, popNames = pops)

# remove NAs from DataFrame
FstDataFrame_noNAs <- na.omit(FstDataFrame)

# identify outliers
outflank_output <- OutFLANK(FstDataFrame_noNAs, 
                            LeftTrimFraction=0.05, 
                            RightTrimFraction=0.05, 
                            Hmin=0.1, 
                            NumberOfSamples=num_pops, 
                            qthreshold=0.05)

# Write output to a text file
outflank_results_df <- outflank_output$results
outlier_df <- filter(outflank_results_df, OutlierFlag == "TRUE")
write.table(outlier_df, file = paste(corename,"_outflank_outliers.txt", ""), sep = '\t', quote = FALSE, row.names = FALSE)

# Plotting results

# plots the actual (yellow) and theoretical (smoothed blue curve) distribution of Fst
OutFLANKResultsPlotter(outflank_output)

# plots Fst against expected Heterozygosity
plot(outflank_output$results$FST, outflank_output$results$He, xlab = "Per Locus FST", ylab = "Per Locus He")

# plots Fst against expected Heterozygosity for outliers
plot(outlier_df$FST, outlier_df$He, xlab = "Per Locus FST", ylab = "Per Locus He")


OF_res_df_logq <- outflank_results_df %>%
  rowwise() %>%
  mutate(logq = log10(qvalues))

# make new plot
plot_log10q_fst <- ggplot()+
  geom_point(data = OF_res_df_logq, aes(x = logq, 
                                 y = FST, color=OutlierFlag), 
             alpha = 0.4)+
  ylab(expression(italic(F[ST]))) +
  scale_color_manual(values=c("black", "red")) +
  xlab("log10 q-value") +
  
  theme_classic() +
  scale_x_reverse() +
  #xlim(0,-6) +
  theme(legend.position="none") #+

plot_log10q_fst

