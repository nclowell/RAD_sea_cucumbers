########################################################################################

# Gene annotation
# 20201218 NL

# get complete list of putatively adaptive loci
# - identified with FST outliers
# - identified with BayEnv2
# - identified with RDA

# make into fasta, and blast against uniprot

# for final df, have column (T/F) for whether this locus was detected
# with these three methods and whether had BLAST results

# provide some summary stats

########################################################################################
### set workding directory and name for naming output files
# set wd
setwd("C:/Users/Natalie Lowell/SHARED_FOLDER/gene_annotation/gene_annotation_202012/")
corename <- "PCtest"

########################################################################################
### input files
BE2_results_path <- "PC_BE2_noEI_likeX_putadap_results.txt"
GO_slim_db_path <- "E:/GO_slim_sorted_db.txt"
blast_results_path <- "PC_putadapt_all3M_1k_regs-uniprot-blastx.tab"


########################################################################################
### import modules
library(tidyverse)
library(hash)

########################################################################################
### functions

# add the ".1" back into SNP names(bad nomenclature decision on my part led to errors downstream)
add_dot1 <- function(snp) {
  namelist <- strsplit(snp,"_")
  newname <- paste0(namelist[[1]][1],".1_",namelist[[1]][2])
  return(newname)
}

# mark Fst outlier SNPs T / F
is_FstOutlier <- function(Locus){
  if(any(fst_outlier_list == Locus)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# mark bayenv2 SNPs T / F
is_besnp <- function(Locus){
  if(any(be2_snp_list == Locus)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# mark RDA SNPs T / F
is_RDAsnp <- function(Locus){
  if(any(rda_snp_list == Locus)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# separate out biolgoical processes IDs into just Ids sep by space
unfold_GOp <- function(GOp){
  if(grepl("[", GOp, fixed=TRUE) == TRUE){ # if GO term, process; else return Nan
    first_list <- strsplit(GOp,";") # list, with "blabla [GO:01]" " blablabla [GO:02]"
    num_items = length(first_list[[1]]) # 2
    new_cell = ""
    for(i in 1:num_items){
      GO <- strsplit(first_list[[1]][i],"[",fixed=TRUE)[[1]][2]
      GO <- strsplit(GO, "]", fixed=TRUE)[[1]][1]
      new_cell <- paste(new_cell,GO)
    }
    
  } else {
    return(NA)
  } # end of else
  new_cell <- trimws(new_cell, which = c("both"))
  return(new_cell)
} # end of fn

# get slim term from id
getSlim <- function(GO_ID, goslim_dict){
  if(is.na(GO_ID)==TRUE || is.null(GO_ID)==TRUE){
    return(NA)
  } else {
    slim <- toString(goslim_dict[[GO_ID]])
    return(slim)
  }
}

########################################################################################
### BayEnv2 results

be_results <- read_tsv(BE2_results_path) %>%
  filter(BayesFactor >= 10) %>% # retain matches with Bayes Factor of at least 10 ("strong" and "decisive")
  select(-X1) # get rid of extra rownames
  
# write to file
write_csv(be_results,paste0(corename,"_BE2_results_atleastBF10.txt"))

# average number of env predictors correlated to each significant SNP
vars_per_locus <- be_results %>%
  count(Locus, EnvVar) %>%
  count(Locus) 

# get number of significant SNPs, average vars per locus, sd vars per locus, for STRONG snps
num_SNPs_sig <- as.numeric(tally(vars_per_locus)[1])
mean_vars_per_SNP <- as.numeric(summarise(vars_per_locus, mean(n)))
sd_vars_per_snp <- as.numeric(summarise(vars_per_locus, sd(n)))
min_vars_per_snp <- as.numeric(summarise(vars_per_locus, min(n)))
max_vars_per_snp <- as.numeric(summarise(vars_per_locus, max(n)))

# how many unique env vars w correlated SNPs? (all 29 env vars)
unique_vars <- be_results %>%
  distinct(EnvVar) %>%
  tally()
unique_vars

be2_snp_list <- vector(length=length(unique(be_results$Locus))) # needs dot1
for(ii in seq(length(be2_snp_list))){
  be2_snp_list[ii] <- add_dot1(unique(be_results$Locus)[ii])
}                   

########################################################################################
### RDA SNPs

# similar df to BE2, with SNP and RDA axis ?
rda_snps <- read_csv("RDA_SNPs_for_annotation.csv") %>%
  rowwise() %>%
  mutate(rda_by_axis = paste0(RDA,"_",axis))

# filter by sig RDAs (2), first axis for each
sig_RDAs <- c("rda_ss_ns", "rda_cv_temp_ns")
rda_snps <- rda_snps %>%
  filter(RDA %in% sig_RDAs) %>%
  filter(axis == 1) # only first axis was significant in significant full model RDAs

rda_snp_list <- unique(rda_snps$snp)

########################################################################################
### FST outlier SNPs

# read in Fst outliers
fst_outliers <- read_csv("PC_Fst_outliers_all_48L.txt", col_names=FALSE)

fst_outliers <- fst_outliers %>%
  rowwise() %>%
  mutate(LocusDotOne = add_dot1(X1))

fst_outlier_list <- fst_outliers$LocusDotOne

########################################################################################
### Get all SNPs across methods (outliers, BE2, RDA)

fst_outlier_list
rda_snp_list
be2_snp_list 

across_methods_unique <- unique(c(fst_outlier_list,rda_snp_list,be2_snp_list ))
write_tsv(as.data.frame(across_methods_unique), 
          paste0(corename, "_putadapt_allMethods_loci.csv"),col_names = FALSE)
length(across_methods_unique)

########################################################################################
### Make GO slim dictionary to translate GO identifier to generic GO slim term for biological processes

# read in goslim generic key
goslim <- read_tsv(GO_slim_db_path, col_names=c("GO_ID", "GO_term","GO_slim_term","Category"))
goslim[1:5,]

# retain only biological processes (Category == P), GO_ID, and high level (GO_slim) terms
goslim_p <- goslim %>%
  filter(Category=="P") %>%
  select(GO_ID,GO_slim_term)
goslim_p[1:5,]
num_entries <- nrow(goslim_p) # get number of rows / entries in goslim database, after filtering for only biological processes

# make dictionary, with keys = GO identifiers and values = GO slim term
goslim_dict <- hash() 
for(entry in 1:num_entries){
  GO_id <- toString(goslim_p[entry,1])
  slim_term <- toString(goslim_p[entry,2])
  goslim_dict[[GO_id]] <- slim_term
}

# test
goslim_dict[["GO:0000001"]]

########################################################################################
### Read in Uniprot GO info (gene ontology per uniprot entry)

better_col_names <- c("Entry","EntryName","Stats","ProteinNames","GeneNames","Organism",
                      "Length","GO_P","GO_C","GO_F","GO_ID")
uniprot_go <- read_tsv("uniprot-reviewed_yes.tab", col_names = better_col_names)
uniprot_go[1:5,]
uniprot_go <- uniprot_go[-c(1),]  # get rid of duplicated old header; dunno how to do this the right way and don't have time rn
uniprot_go[1:5,]

# keep only biological processes GO
uniprot_go_p <- select(uniprot_go,-c(GO_C, GO_F, GO_ID))
uniprot_go_p[1:5,]

########################################################################################
### Read in blastx results (our matches to the uniprot blast database)

# add informative format 6 blastx results header to results files that are currently headerless
fmt6_header = c("qseqid", "sseqid", "pident", "length", 
               "mismatch", "gapopen", "qstart", "qend",
               "sstart", "send", "evalue", "bitscore")
blast <- read_tsv(blast_results_path, col_names = fmt6_header)


########################################################################################
### Merge gene ontology information & get GO slim terms

# in order to merge, extract Entry from sseqid first
blast_wEntry <- blast %>%
  rowwise() %>%
  mutate(Entry=strsplit(sseqid,"|",fixed=TRUE)[[1]][2])
blast_wGOp <- left_join(blast_wEntry, uniprot_go_p, by = "Entry")

########################################################################################
### Merge blastx results with uniprot GO_p terms

# pull out GO IDs for biological processes and create row for each unique GO ID
blast_wGOp_unfold <- blast_wGOp  %>%
  rowwise() %>%
  mutate(GO_ID_p=unfold_GOp(GO_P)) %>%
  separate_rows(GO_ID_p, sep=" ")

# add slim column, with broader GO categories
blast_wGOp_unfold_wslim <- blast_wGOp_unfold %>%
  rowwise() %>%
  mutate(GO_slim = getSlim(GO_ID=GO_ID_p, goslim_dict=goslim_dict))


########################################################################################
### Add categories for whether FST outlier, BE2 SNP, RDA SNP

# mark SNPs based on method

blast_wGOp_unfold_wslim_1M <- blast_wGOp_unfold_wslim %>%
  rowwise() %>%
  mutate(IS_FstOutSNP = is_FstOutlier(qseqid))

blast_wGOp_unfold_wslim_2M <- blast_wGOp_unfold_wslim_1M %>%
  rowwise() %>%
  mutate(IS_BE2SNP = is_besnp(qseqid))

blast_wGOp_unfold_wslim_3M <- blast_wGOp_unfold_wslim_2M %>%
  rowwise() %>%
  mutate(IS_RDASNP = is_RDAsnp(qseqid))

########################################################################################
### Merge env predictor info with blastx GO slim data, from BE2

be_results_dot1 <- be_results %>%
  rowwise() %>%
  mutate(Locus = add_dot1(Locus))

blast_wGOp_unfold_wslim_3M_pred <- full_join(be_results_dot1, blast_wGOp_unfold_wslim_3M,by=c("Locus" = "qseqid"))

# add flag for blast hit
blast_wGOp_unfold_wslim_3M_pred_bl <- blast_wGOp_unfold_wslim_3M_pred %>%
  rowwise() %>%
  mutate(BLASThit = !is.na(sseqid))

########################################################################################
### BayEnv2 questions

# How many significant SNPs per env predictor var?
BE2_snps_per_var <- be_results %>%
  group_by(EnvVar) %>%
  tally() %>%
  arrange(desc(n))
write_csv(BE2_snps_per_var, paste0(corename, "_BE2_snps_per_var.csv"))

# What are the common bio processes per env predictor var?
BE2_env_GOslimP_counts <- blast_wGOp_unfold_wslim_3M_pred_bl %>%
  filter(IS_BE2SNP == TRUE) %>%
  filter(BLASThit == TRUE) %>%
  group_by(EnvVar, GO_slim) %>%
  tally() %>%
  arrange(EnvVar, desc(n)) %>%
  drop_na(GO_slim) %>%
  filter(GO_slim != "")
write_csv(BE2_env_GOslimP_counts,
          paste0(corename, "_BE2_env_GOslimP_counts.csv"))

# What is the proportion of each bio process per env predictor var?
BE2_env_GOslimP_props <- blast_wGOp_unfold_wslim_3M_pred_bl %>%
  filter(IS_BE2SNP) %>% # must be a BayEnv2 SNP
  filter(BLASThit) %>% # must have GO annotation
  group_by(EnvVar, GO_slim) %>%
  tally() %>%
  drop_na(GO_slim) %>% 
  filter(GO_slim != "") %>%
  filter(GO_slim != "other metabolic processes") %>%
  filter(GO_slim != "other biological processes") %>%
  arrange(EnvVar, desc(n)) %>%
  group_by(EnvVar) %>%
  #summarise(n = n()) %>%
  mutate(freq = n / sum(n), by = EnvVar) %>%
  select(-by)
write_csv(BE2_env_GOslimP_props,
          paste0(corename, "_BE2_env_GOslimP_props.csv"))

# How many SNPs have GO slim annotations per env var?
num_BEsnp_withGO_by_var <- blast_wGOp_unfold_wslim_3M_pred_bl %>%
  filter(IS_BE2SNP) %>% # must be a BayEnv2 SNP
  filter(BLASThit) %>% # must have GO annotation
  group_by(EnvVar, Locus) %>%
  tally() %>%
  count(EnvVar)
write_csv(num_BEsnp_withGO_by_var,
          "_num_BEsnp_withGO_by_var_20200311.csv")

# rank biological processes per annotated snp per env var
goslims_X_BEsnp_X_var <- blast_wGOp_unfold_wslim_3M_pred_bl %>%
  filter(IS_BE2SNP) %>% # must be a BayEnv2 SNP
  filter(BLASThit) %>% # must have GO annotation
  group_by(EnvVar, GO_slim) %>%
  tally() %>%
  drop_na(GO_slim) %>% 
  filter(GO_slim != "") %>%
  filter(GO_slim != "other metabolic processes") %>%
  filter(GO_slim != "other biological processes") %>%
  arrange(desc(n))
write_csv(goslims_X_BEsnp_X_var, paste0(corename,"_goslims_X_BEsnp_X_var.csv"))


# proportions instead, by env var
rank_GOslims_by_var_props <- blast_wGOp_unfold_wslim_3M_pred_bl %>%
  filter(IS_BE2SNP) %>% # must be a BayEnv2 SNP
  filter(BLASThit) %>% # must have GO annotation
  group_by(EnvVar, GO_slim) %>%
  tally() %>%
  drop_na(GO_slim) %>% 
  filter(GO_slim != "") %>%
  filter(GO_slim != "other metabolic processes") %>%
  filter(GO_slim != "other biological processes") %>%
  arrange(EnvVar, desc(n)) %>%
  group_by(EnvVar) %>%
  #summarise(n = n()) %>%
  mutate(freq = n / sum(n), by = EnvVar) 
write_csv(rank_GOslims_by_var_props, paste0(corename,"_rank_GOslims_by_var_props.csv"))

##############################################################
### Other questions

# How many unique snps from bayenv had blast hits? 15
num_snps_besnp_wBLAST <- blast_wGOp_unfold_wslim_3M_pred_bl %>%
  filter(IS_BE2SNP) %>% # must be a BayEnv2 SNP
  filter(BLASThit) %>% # must have GO annotation
  group_by(Locus) %>%
  tally() %>%
  count()

# How many unique snps from fst outlers had blast hits? 6
num_snps_fstoutlier_wBLAST <- blast_wGOp_unfold_wslim_3M_pred_bl %>%
  filter(IS_FstOutSNP) %>% # must be a BayEnv2 SNP
  filter(BLASThit) %>% # must have GO annotation
  group_by(Locus) %>%
  tally() %>%
  count()

# How many unique snps from RDA outlers had blast hits? 6
num_snps_RDA_wBLAST <- blast_wGOp_unfold_wslim_3M_pred_bl %>%
  filter(IS_RDASNP) %>% # must be a BayEnv2 SNP
  filter(BLASThit) %>% # must have GO annotation
  group_by(Locus) %>%
  tally() %>%
  count()


##############################################################
### FST Outliers

# Get biological processes for SNPs that are only FST outliers
FSTonly_goslims <- blast_wGOp_unfold_wslim_3M_pred_bl %>%
  filter(IS_FstOutSNP) %>%
  filter(!IS_BE2SNP) %>%
  filter(!IS_RDASNP) %>%
  filter(BLASThit) %>%
  group_by(GO_slim) %>%
  tally() %>%
  arrange(GO_slim, desc(n)) %>% # not working?
  filter(GO_slim != "")
write_csv(FSTonly_goslims,
          paste0(corename, "_FSTonly_goslims.csv"))


##############################################################
### Curious as to whether results are based on random selection
# from abundance of different goslim p terms

# What is the frequency of GO Slim terms in the database?
goslim_p_counts <- goslim_p %>%
  group_by(GO_slim_term) %>%
  tally() %>%
  arrange(desc(n))

write_csv(goslim_p_counts,
          "goslim_database_abundaces.csv")


# Get biological processes for SNPs that are only FST outliers
FSTonly_goslims <- blast_wGOp_unfold_wslim_3M_pred_bl %>%
  filter(IS_FstOutSNP) %>%
  filter(!IS_BE2SNP) %>%
  filter(!IS_RDASNP) %>%
  filter(BLASThit) %>%
  group_by(GO_slim) %>%
  tally() %>%
  arrange(GO_slim, desc(n)) %>% # not working?
  filter(GO_slim != "")
write_csv(FSTonly_goslims,
          paste0(corename, "_FSTonly_goslims.csv"))





