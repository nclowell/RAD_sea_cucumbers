################################################################################
### Get He, Ho, afs, and Frequency polymorphic loci per population from a genepop file
# 20200327 Natalie Lowell nclowell@uw.edu
# 
# Assumptions:
# - Individual names are population underscore individual identifier, like
# "AukeBay_25"
# - Using Version of strataG 2.4.905 (objects changed a lot through versions)

################################################################################

# import libraries
library(adegenet)
library(strataG)
library(tidyverse)

################################################################################
### Edit these things first, then run code in next section

# set working directory; this is where output files will be stored
setwd("E:/dDocent_for_mox/parcal_wd/parcal_mox001/analyses_CG_PC/prop_polym")

# path to gp file
gpfilepath <- "E:/dDocent_for_mox/parcal_wd/parcal_mox001/parcal_mox001_qc/parcal_mox001_qc.gen"

# number of digits coding single genotype in genepop file (often 2 or 3)
my_ncode <- 3

# set "core" name to be assigned to output files
corename <- "PC_all_20200330"

################################################################################

# read in genepop file
genind <-read.genepop(gpfilepath, ncode = 3)
gtype <- genind2gtypes(genind)
gtype_data_id <- gtype@data$id
inds <- unique(gtype_data_id)
gtype_data_locus <- gtype@data$locus
locus_names <- unique(gtype_data_locus)
numloci <- length(locus_names)

# write function to get population name from individual name
get_pop <- function(ind){
  ind <- toString(ind)
  pop <- strsplit(ind,"_")[[1]][1]
  return(pop)
}

# get population for each individual, assuming inds are named
# pop underscore ind identifier, with no underscores anywhere else;
# assign to strata
pop_strata <- vector(length=length(gtype_data_id))
for(j in seq(1,length(gtype_data_id))){
  pop_strata[j] <- get_pop(gtype_data_id[j])
}
popnames <- unique(pop_strata)
numpops <- length(popnames)
gtype@data$stratum <- pop_strata

### Allele frequencies per population

# function for getting allele frequencies by
make_afs_df <- function(gtype){
  afs_df <- data.frame("LocusName" = character(), "Pop" = character(),
                       "AlleleName" = character(), "Frequency" = numeric(),
                       stringsAsFactors = FALSE)
  afs <- alleleFreqs(gtype, by.strata = TRUE, type = "prop")
  rowcount = 1
  for(locus_index in seq(1,numloci)){
    for(popnum in seq(1, numpops)){
      for(allele in seq(1,nrow(afs[[locus_index]]))){
        prop_allele <- afs[[locus_index]][allele,popnum]
        alleleName <- dimnames(afs[[locus_index]])[[1]][allele]
        afs_df[rowcount,] <- c(locus_names[locus_index], popnames[popnum], alleleName, prop_allele)
        rowcount <- rowcount + 1
      }
    }
  }
  return(afs_df)
}
afs_by_pop <- make_afs_df(gtype)
afs_tidy <- as_tibble(afs_by_pop)
afs_tidy$Frequency <- as.numeric(as.character(afs_tidy$Frequency))
class(afs_tidy$Frequency)

is_polym <- function(af){
  if(is.na(af)){
    return(NA)
  } else {
    if(af == 1 | af == 0){
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
}

# add column for whether locus polymoprhic per population
afs_tidy_plm <- afs_tidy %>%
  rowwise() %>%
  mutate(polym = is_polym(Frequency))

# (locus) counts T / F polymorphic per population
counts_polym <- afs_tidy_plm %>%
  group_by(Pop) %>%
  count(polym)

# total observations
totals <- counts_polym %>%
  ungroup() %>%
  group_by(Pop) %>%
  summarise(total = sum(n))

counts_totals_polym <- full_join(counts_polym, totals, by = "Pop") %>%
  filter(polym == TRUE) %>%
  rowwise() %>%
  mutate(prop_polym = n / total) %>%
  select(-polym, -n, -total)

### Heterozygosity, per population

# table
He <- heterozygosity(gtype, by.strata = TRUE, type = "expected")
He_tidy <- as_tibble(He)
He_means <- He_tidy %>%
  group_by(stratum) %>%
  summarise(mean(exptd.het))
Ho <- heterozygosity(gtype, by.strata = TRUE, type = "observed")
Ho_tidy <- as_tibble(Ho)
Ho_means <- Ho_tidy %>%
  group_by(stratum) %>%
  summarise(mean(obsvd.het))

results_H <- full_join(He_means, Ho_means, by = c("stratum"))
results_H_pP <- full_join(results_H, counts_totals_polym, by = c("stratum" = "Pop")) 
names(results_H_pP) <- c("Pop", "Mean He", "Mean Ho", "Frequency polymorphic")
write_csv(results_H_pP, paste0(corename, "_He_Ho_propPolym_results.csv"))

# plots
Ho_plot <- ggplot(data=Ho_tidy, aes(x=obsvd.het)) +
  geom_histogram() +
  theme_classic() +
  labs(x="Observed heterozygosity (Ho)", y= "Count") + 
  facet_wrap(~stratum)
Ho_plot_name <- paste0(corename, "_Ho_hists.png")
ggsave(Ho_plot_name, Ho_plot, height = 6, width = 4)

He_plot <- ggplot(data=He_tidy, aes(x=exptd.het)) +
  geom_histogram() +
  theme_classic() +
  labs(x="Expected heterozygosity (He)", y= "Count") + 
  facet_wrap(~stratum)
He_plot_name <- paste0(corename, "_He_hists.png")
ggsave(He_plot_name, He_plot, height = 6, width = 4)
  