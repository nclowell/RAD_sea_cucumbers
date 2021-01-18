### IBD nuetral & all loci PARCAL only
# 20200803 NL
# and looking for corr with linear km and log km

#################################################################################
# get packages into R
library(geosphere)
library(dplyr)
library(viridis)
library(ggplot2)
library(vegan)
library(tidyverse)

setwd("E:/dDocent_for_mox/parcal_wd/parcal_mox001/IBD/")

# read in file with lat & long, pairwise FST, etc.
parcal_data <- read.csv("PC_all_neu_adap_FST_inwater_20201016.csv")

#################################################################################
# mantel takes two arrays
# - (1) pairwise distance matrix, as the crow flies and also in water
# - (2) pairwise FST matrix, here one for all loci and also one for neutral loci

# save the character values of the population names
pop_name <- with(parcal_data, sort(unique(c(as.character(pop1),
                                            as.character(pop2)))))

#create some  empty 2-D  arrays to hold the pairwise data
dist_array<- array(data = 0, dim = c(length(pop_name), length(pop_name)), 
                   dimnames = list(pop_name, pop_name))

fst_array <- array(data = 0, dim = c(length(pop_name), length(pop_name)), 
                   dimnames = list(pop_name, pop_name))
neu_fst_array <- array(data = 0, dim = c(length(pop_name), length(pop_name)), 
                   dimnames = list(pop_name, pop_name))
adap_fst_array <- array(data = 0, dim = c(length(pop_name), length(pop_name)), 
                       dimnames = list(pop_name, pop_name))
# crow_dist_array <- array(data = 0, dim = c(length(pop_name), length(pop_name)), 
#                          dimnames = list(pop_name, pop_name))
# logcrow_array <- array(data = 0, dim = c(length(pop_name), length(pop_name)), 
#                          dimnames = list(pop_name, pop_name))
# log_dist_array<- array(data = 0, dim = c(length(pop_name), length(pop_name)), 
#                    dimnames = list(pop_name, pop_name))

# save some vectors of the positions of first matches of the first argument to the second
i <- match(parcal_data$pop1, pop_name)
j <- match(parcal_data$pop2, pop_name)

# linearize Fst, add column to dataframe
parcal_data <- parcal_data %>% # ALL LOCI
  rowwise() %>%
  mutate(linearized_fst = (fst/(1-fst)))

# linearize Fst, add column to dataframe
parcal_data <- parcal_data %>% # NEUTRAL LOCI
  rowwise() %>%
  mutate(linearized_neu_fst = (neu_fst/(1-neu_fst)))

parcal_data <- parcal_data %>% # NEUTRAL LOCI
  rowwise() %>%
  mutate(linearized_adap_fst = (adap_fst/(1-adap_fst)))

# add as crow flies
# parcal_data <- parcal_data %>% 
#   rowwise() %>%
#   mutate(km2 = (distVincentyEllipsoid(c(pop1long, pop1lat), c(pop2long, pop2lat))/1000)) %>%
#   rowwise() %>%
#   mutate(log_crow_km = log(km2))

# make inwater distance as numeric
parcal_data$inwater_dist_km = as.numeric(parcal_data$inwater_dist_km)
# parcal_data$log_km = as.numeric(parcal_data$log_km)
# parcal_data$log_crow_km = as.numeric(parcal_data$log_crow_km)

# populate the empty arrays with data saved in the vectors
dist_array[cbind(i,j)] <- dist_array[cbind(j,i)] <- parcal_data$inwater_dist_km
#log_dist_array[cbind(i,j)] <- log_dist_array[cbind(j,i)] <- parcal_data$log_km
#crow_dist_array[cbind(i,j)] <- crow_dist_array[cbind(j,i)] <- parcal_data$km2
#logcrow_array[cbind(i,j)] <- logcrow_array[cbind(j,i)] <- parcal_data$log_crow_km
fst_array[cbind(i,j)] <- fst_array[cbind(j,i)] <- parcal_data$linearized_fst
neu_fst_array[cbind(i,j)] <- neu_fst_array[cbind(j,i)] <- parcal_data$linearized_neu_fst
adap_fst_array[cbind(i,j)] <- adap_fst_array[cbind(j,i)] <- parcal_data$linearized_adap_fst

write.csv(dist_array, "inwater_dist_array_20201016.csv")
#write.csv(log_dist_array, "log_dist_array_20200803.csv")
#write.csv(logcrow_array, "log_crow_array_20200803.csv")

#write.csv(crow_dist_array, "crow_dist_array_20200630.csv")
write.csv(fst_array, "all_loci_fst_array_20201016.csv")
write.csv(neu_fst_array, "neu_loci_fst_array_20201016.csv")
write.csv(adap_fst_array, "adpa_loci_fst_array_20201016.csv")

######### run mantel test, all loci vs. inwater
mantel(dist_array, fst_array, method="pearson", permutations=10000)

# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = dist_array, ydis = fst_array, method = "pearson",      permutations = 10000) 
# 
# Mantel statistic r: 0.7581 
#       Significance: 0.00019998 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.219 0.317 0.410 0.504 
# Permutation: free
# Number of permutations: 10000

# plot
ggplot(data = parcal_data, aes(x=inwater_dist_km, y=linearized_fst)) + #specify dataframe
  geom_point( size = 3, alpha = 0.9) +
  ylab(expression(italic(F[ST]/(1-F[ST])))) +                           #set labels for the axes and title
  xlab("In-water distance (km)")  +
  geom_smooth(method = "lm") + # add regression line
  theme_classic() +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size=16))
ggsave(filename = "IBD_parcal_allloci_inwater_20201016.png", dpi=300)

lm_results <- lm(parcal_data$linearized_fst~parcal_data$inwater_dist_km)
lm_results$coefficients
# (Intercept) parcal_data$inwater_dist_km 
# 3.492701e-03                3.184471e-06 

######### run mantel test,  neutral, inwater
mantel(dist_array, neu_fst_array, method="pearson", permutations=10000)
# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = dist_array, ydis = neu_fst_array, method = "pearson",      permutations = 10000) 
# 
# Mantel statistic r: 0.4891 
# Significance: 0.0058994 
# 
# Upper quantiles of permutations (null model):
# 90%   95% 97.5%   99% 
# 0.246 0.323 0.382 0.453 
# Permutation: free
# Number of permutations: 10000

# plot
ggplot(data = parcal_data, aes(x=inwater_dist_km, y=linearized_neu_fst)) + #specify dataframe
  geom_point( size = 3, alpha = 0.9) +
  ylab(expression(italic(F[ST]/(1-F[ST])))) +                           #set labels for the axes and title
  xlab("In-water distance (km)")  +
  geom_smooth(method = "lm") + # add regression line
  theme_classic() +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size=16))
ggsave(filename = "IBD_parcal_neuloci_inwater_20201016.png", dpi=300)

lm_results <- lm(parcal_data$linearized_neu_fst~parcal_data$inwater_dist_km)
lm_results$coefficients
# (Intercept) parcal_data$inwater_dist_km 
# 1.786348e-03                1.079283e-06 


######### run mantel test,  adap, inwater
mantel(dist_array, adap_fst_array, method="pearson", permutations=10000)
# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = dist_array, ydis = adap_fst_array, method = "pearson",      permutations = 10000) 
# 
# Mantel statistic r: 0.7855 
#       Significance: 9.999e-05 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.220 0.299 0.383 0.473 
# Permutation: free
# Number of permutations: 10000

# plot
ggplot(data = parcal_data, aes(x=inwater_dist_km, y=linearized_adap_fst)) + #specify dataframe
  geom_point( size = 3, alpha = 0.9) +
  ylab(expression(italic(F[ST]/(1-F[ST])))) +                           #set labels for the axes and title
  xlab("In-water distance (km)")  +
  geom_smooth(method = "lm") + # add regression line
  theme_classic() +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size=16))
ggsave(filename = "IBD_parcal_adaploci_inwater_20201016.png", dpi=300)

lm_results <- lm(parcal_data$linearized_adap_fst~parcal_data$inwater_dist_km)
lm_results$coefficients
# (Intercept) parcal_data$inwater_dist_km 
# 1.898837e-02                2.196213e-05 


### layer all IBD plots??

fix <- function(type){
  if(type=="linearized_fst"){
    return("All")}
  if(type=="linearized_neu_fst"){
    return("Neutral")}
  if(type=="linearized_adap_fst"){
    return("Adaptive") }
}

tidyfst <- parcal_data %>%
  select(pop1,pop2,inwater_dist_km,linearized_fst, linearized_adap_fst, linearized_neu_fst) %>%
  gather(4:6, key="which_loci", value="linearized_fst") %>%
  rowwise() %>%
  mutate(Type=fix(which_loci)) %>%
  select(-which_loci)

ggplot(data = tidyfst, aes(x=inwater_dist_km, y=linearized_fst, color = Type)) + #specify dataframe
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm") + # add regression line
  scale_colour_manual(values=c("red", "black", "blue")) +
  ylab(expression(italic(F[ST]/(1-F[ST])))) + #set labels for the axes and title
  xlab("In-water distance (km)")  +
  theme_classic() +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size=16)) +
  theme(legend.position="bottom") +
  scale_fill_discrete(labels = c("Adaptive", "All", "Neutral")) +
  theme(legend.title=element_blank())

#ggsave(filename = "IBD_parcal_adaploci_inwater_20201016.png", dpi=300)


