####################################################################
### AMOVA in R
# 20201014 Natalie Lowell

####################################################################

# import libraries
library(poppr)
library(pegas)

####################################################################

# set working directory 
setwd("E:/dD_cragig_mox_20191205/filtering/analyses/amova/pegas")

####################################################################

# name & read in genepop files, all loci, neutral, and putatively adaptive
PC_gp_file_all = "parcal_mox001_qc.gen"
PC_gp_file_neu = "PC_likeX_neu1877.gen"
PC_gp_file_adap = "PC_likeX_putadapt198.gen"

####################################################################

PC_data_all <-read.genepop(PC_gp_file_all, ncode=3)
PC_data_neu <-read.genepop(PC_gp_file_neu, ncode=3)
PC_data_adap <-read.genepop(PC_gp_file_adap, ncode=3)

# read in strata file
PC_strata<- read.csv("PC_strata_for_amova.csv")
head(PC_strata)

# assign hierarchical levels from strata file to genind object
strata(PC_data_all) <- PC_strata
strata(PC_data_neu) <- PC_strata
strata(PC_data_adap) <- PC_strata

# change the format of your data from a genind to a genclone object
PC_all_gc <- as.genclone(PC_data_all)
PC_neu_gc <- as.genclone(PC_data_neu)
PC_adap_gc <- as.genclone(PC_data_adap)

amova_PC_site_all_nowithin <- poppr.amova(PC_all_gc, ~site, within = FALSE, method = "pegas", nperm = 1000)
amova_PC_site_all_nowithin

amova_PC_site_neu_nowithin <- poppr.amova(PC_neu_gc, ~site, within = FALSE, method = "pegas", nperm = 1000)
amova_PC_site_neu_nowithin

amova_PC_site_adap_nowithin <- poppr.amova(PC_adap_gc, ~site, within = FALSE, method = "pegas", nperm = 1000)
amova_PC_site_adap_nowithin

# by state
amova_PC_state_all_nowithin <- poppr.amova(PC_all_gc, ~state, within = FALSE, method = "pegas", nperm = 1000)
amova_PC_state_all_nowithin

amova_PC_state_neu_nowithin <- poppr.amova(PC_neu_gc, ~state, within = FALSE, method = "pegas", nperm = 1000)
amova_PC_state_neu_nowithin

amova_PC_state_adap_nowithin <- poppr.amova(PC_adap_gc, ~state, within = FALSE, method = "pegas", nperm = 1000)
amova_PC_state_adap_nowithin

# by PS
amova_PC_PS_all_nowithin <- poppr.amova(PC_all_gc, ~PS, within = FALSE, method = "pegas", nperm = 1000)
amova_PC_PS_all_nowithin

amova_PC_PS_neu_nowithin <- poppr.amova(PC_neu_gc, ~PS, within = FALSE, method = "pegas", nperm = 1000)
amova_PC_PS_neu_nowithin

amova_PC_PS_adap_nowithin <- poppr.amova(PC_adap_gc, ~PS, within = FALSE, method = "pegas", nperm = 1000)
amova_PC_PS_adap_nowithin







