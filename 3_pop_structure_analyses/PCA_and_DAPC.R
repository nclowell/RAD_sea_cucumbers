### PCA & DAPC for all loci, putatively neutral loci, and putatively adaptive loci
# 20200310 NL

####

library(adegenet)
library (ggplot2)
library(ggpubr)

####

setwd("E:/dDocent_for_mox/parcal_wd/parcal_mox001/analyses_CG_PC/bayenv2/BayEnv2_round2_likeX/PC_CG_neu_adap/")

# Read in genepop files of neutral loci loci only and outliers only
PC_gpfile_all <-read.genepop("E:/dDocent_for_mox/parcal_wd/parcal_mox001/parcal_mox001_qc/parcal_mox001_qc.gen", ncode = 3)
PC_gpfile_neu <-read.genepop("PC_likeX_neu1877.gen", ncode = 3)
PC_gpfile_out <-read.genepop("PC_likex_putadapt198.gen", ncode = 3)

# Replace missing data information with the mean
PC_gpfile_all_scaled <- scaleGen(PC_gpfile_all, NA.method="mean")
PC_gpfile_neu_scaled <- scaleGen(PC_gpfile_neu, NA.method="mean")
PC_gpfile_out_scaled <- scaleGen(PC_gpfile_out, NA.method="mean")


PC_ordered_pops <- c("Chiniak Bay, AK",
                     "Yakutat Bay, AK",  
                     "Auke Bay, AK",   
                     "Bella Bella, BC",
                     "Sekiu, WA",
                     "James Island, WA",
                     "Keyport, WA",
                     "Eld Inlet, WA",
                     "Charleston, OR")

PC_site_colors <- c("darkorchid4",
                    "darkblue",
                    "cadetblue",
                    "aquamarine2",
                    "green4",
                    "greenyellow",
                    "gold",
                    "chocolate1",
                    "red4")
  
PC_ordered_states <- c("Alaska",
                       "BritishColumbia",
                       "Washington",
                       "Oregon")

PC_state_colors <- c("tomato",
                    "goldenrod2",
                    "dodgerblue4",
                    "mediumorchid3")

ordered_SS <- c("Inside", "Outside")
SS_colors = c("dodgerblue4", "orange3")

ordered_NS <- c("North", "South")
NS_colors = c("darkorchid4", "goldenrod2")


##############################################################################
### PCAs

# Conduct PCAs
PC_pca_gpfile_all <- dudi.pca(PC_gpfile_all_scaled,cent=TRUE,scale=FALSE,scannf=FALSE,nf=3)
PC_pca_gpfile_neu <- dudi.pca(PC_gpfile_neu_scaled,cent=TRUE,scale=FALSE,scannf=FALSE,nf=3)
PC_pca_gpfile_out <- dudi.pca(PC_gpfile_out_scaled,cent=TRUE,scale=FALSE,scannf=FALSE,nf=3)

# Save the plotting coordinates as a seprate dataframe, so you can customize plot
PC_gpfile_all_df <- PC_pca_gpfile_all$li
PC_gpfile_neu_df <- PC_pca_gpfile_neu$li  
PC_gpfile_out_df <- PC_pca_gpfile_out$li  

PC_popnames_gpfile <- c(rep("Chiniak Bay, AK", 19), 
                     rep("Yakutat Bay, AK", 39),
                     rep("Auke Bay, AK",10),
                     rep("Bella Bella, BC", 32), 
                     rep("Sekiu, WA", 47),
                     rep("James Island, WA",50),
                     rep("Keyport, WA",42),
                     rep("Eld Inlet, WA", 48), 
                     rep("Charleston, OR", 38))
PC_gpfile_all_df$Population = PC_popnames_gpfile # add the names to the df
PC_gpfile_neu_df$Population = PC_popnames_gpfile
PC_gpfile_out_df$Population = PC_popnames_gpfile 

# re-order the levels in the order of appearance in the data.frame
PC_gpfile_all_df$Population <- factor(PC_gpfile_all_df$Population, levels = PC_ordered_pops)
PC_gpfile_neu_df$Population <- factor(PC_gpfile_neu_df$Population, levels = PC_ordered_pops)
PC_gpfile_out_df$Population <- factor(PC_gpfile_out_df$Population, levels = PC_ordered_pops)

PC_statenames_gpfile <-  c(rep("Alaska", 19), 
                        rep("Alaska", 39),
                        rep("Alaska",10),
                        rep("BritishColumbia", 32), 
                        rep("Washington", 47),
                        rep("Washington",50),
                        rep("Washington",42),
                        rep("Washington", 48), 
                        rep("Oregon", 38))
PC_gpfile_all_df$States = PC_statenames_gpfile # add the names to the df
PC_gpfile_neu_df$States = PC_statenames_gpfile 
PC_gpfile_out_df$States = PC_statenames_gpfile 

# re-order the levels in the order of appearance in the data.frame
PC_gpfile_all_df$States <- factor(PC_gpfile_all_df$States, levels = PC_ordered_states)
PC_gpfile_neu_df$States <- factor(PC_gpfile_neu_df$States, levels = PC_ordered_states)
PC_gpfile_out_df$States <- factor(PC_gpfile_out_df$States, levels = PC_ordered_states)
PC_reSS_gpfile <- c(rep("Outside", 19), 
                 rep("Outside", 39),
                 rep("Outside",10),
                 rep("Outside", 32), 
                 rep("Inside", 47),
                 rep("Inside",50),
                 rep("Inside",42),
                 rep("Inside", 48), 
                 rep("Outside", 38))
PC_gpfile_all_df$SalishSea = PC_reSS_gpfile # add the names to the df
PC_gpfile_neu_df$SalishSea = PC_reSS_gpfile 
PC_gpfile_out_df$SalishSea = PC_reSS_gpfile

# re-order the levels in the order of appearance in the data.frame
PC_gpfile_all_df$SalishSea <- factor(PC_gpfile_all_df$SalishSea, levels = ordered_SS)
PC_gpfile_neu_df$SalishSea <- factor(PC_gpfile_neu_df$SalishSea, levels = ordered_SS)
PC_gpfile_out_df$SalishSea <- factor(PC_gpfile_out_df$SalishSea, levels = ordered_SS)

PC_NSnames_gpfile <-  c(rep("North", 19), 
                           rep("North", 39),
                           rep("North",10),
                           rep("North", 32), 
                           rep("South", 47),
                           rep("South",50),
                           rep("South",42),
                           rep("South", 48), 
                           rep("South", 38))

PC_gpfile_all_df$NS = PC_NSnames_gpfile # add the names to the df
PC_gpfile_neu_df$NS = PC_NSnames_gpfile
PC_gpfile_out_df$NS = PC_NSnames_gpfile

# re-order the levels in the order of appearance in the data.frame
PC_gpfile_all_df$SalishSea <- factor(PC_gpfile_all_df$SalishSea, levels = ordered_NS)
PC_gpfile_neu_df$SalishSea <- factor(PC_gpfile_neu_df$SalishSea, levels = ordered_NS)
PC_gpfile_out_df$SalishSea <- factor(PC_gpfile_out_df$SalishSea, levels = ordered_NS)

# get percet var explained by each axis
PC_eigs_all <- PC_pca_gpfile_all$eig
PC_axis1_var_all <- paste(substr(as.character((PC_eigs_all[1] / sum(PC_eigs_all))*100),1,5),"%",sep="")
PC_axis2_var_all <- paste(substr(as.character((PC_eigs_all[2] / sum(PC_eigs_all))*100),1,5),"%",sep="")
PC_eigs_neu <- PC_pca_gpfile_neu$eig
PC_axis1_var_neu <- paste(substr(as.character((PC_eigs_neu[1] / sum(PC_eigs_neu))*100),1,5),"%",sep="")
PC_axis2_var_neu <- paste(substr(as.character((PC_eigs_neu[2] / sum(PC_eigs_neu))*100),1,5),"%",sep="")
PC_eigs_out <- PC_pca_gpfile_out$eig
PC_axis1_var_out <- paste(substr(as.character((PC_eigs_out[1] / sum(PC_eigs_out))*100),1,5),"%",sep="")
PC_axis2_var_out <- paste(substr(as.character((PC_eigs_out[2] / sum(PC_eigs_out))*100),1,5),"%",sep="")

########################
### Plot PCAs

############
# by site

PC_gpfile_all_df$Axis1 <- PC_gpfile_all_df$Axis1*-1 # rotating across axis to make it easier to read
PC_PCA_all_bysite <- ggplot(data = PC_gpfile_all_df, aes(x= Axis1, y= Axis2)) + 
  geom_point(aes(colour = Population), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values=PC_site_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_axis1_var_all, y=PC_axis2_var_all) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ggtitle("P. californicus") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="right") +
  theme(legend.title=element_blank())
PC_PCA_all_bysite

PC_PCA_neu_bysite <- ggplot(data = PC_gpfile_neu_df, aes(x= Axis1, y= Axis2)) + 
  geom_point(aes(colour = Population), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values=PC_site_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_axis1_var_neu, y=PC_axis2_var_neu) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="right") +
  theme(legend.title=element_blank())
PC_PCA_neu_bysite

PC_PCA_out_bysite <- ggplot(data = PC_gpfile_out_df, aes(x= Axis1, y= Axis2)) + 
  geom_point(aes(colour = Population), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values=PC_site_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_axis1_var_out, y=PC_axis2_var_out) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="right") +
  theme(legend.title=element_blank())
PC_PCA_out_bysite

############
# by state

PC_PCA_all_bystate <- ggplot(data = PC_gpfile_all_df, aes(x= Axis1, y= Axis2)) + 
  geom_point(aes(colour = States), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values=PC_state_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_axis1_var_all, y=PC_axis2_var_all) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="right") +
  theme(legend.title=element_blank())
PC_PCA_all_bystate

PC_PCA_neu_bystate <- ggplot(data = PC_gpfile_neu_df, aes(x= Axis1, y= Axis2)) + 
  geom_point(aes(colour = States), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values=PC_state_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_axis1_var_neu, y=PC_axis2_var_neu) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="right") +
  theme(legend.title=element_blank())
PC_PCA_neu_bystate

PC_PCA_out_bystate <- ggplot(data = PC_gpfile_out_df, aes(x= Axis1, y= Axis2)) + 
  geom_point(aes(colour = States), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values=PC_state_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_axis1_var_out, y=PC_axis2_var_out) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="right") +
  theme(legend.title=element_blank())
PC_PCA_out_bystate

############
# by reSS

PC_PCA_all_bySS <- ggplot(data = PC_gpfile_all_df, aes(x= Axis1, y= Axis2)) + 
  geom_point(aes(colour = SalishSea), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values = SS_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_axis1_var_all, y=PC_axis2_var_all) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="right") +
  theme(legend.title=element_blank())
PC_PCA_all_bySS

PC_PCA_neu_bySS <- ggplot(data = PC_gpfile_neu_df, aes(x= Axis1, y= Axis2)) + 
  geom_point(aes(colour = SalishSea), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values = SS_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_axis1_var_neu, y=PC_axis2_var_neu) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="right") +
  theme(legend.title=element_blank())
PC_PCA_neu_bySS

PC_PCA_out_bySS <- ggplot(data = PC_gpfile_out_df, aes(x= Axis1, y= Axis2)) + 
  geom_point(aes(colour = SalishSea), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values=SS_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_axis1_var_out, y=PC_axis2_var_out) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="right") +
  theme(legend.title=element_blank())
PC_PCA_out_bySS

##############################################################################
### DAPCs

PC_ninds = 325
PC_npcas = PC_ninds - 1
PC_n_pops = 9
PC_ndims = PC_n_pops - 1

PC_dapc_all <- dapc(PC_gpfile_all,PC_gpfile_all$pop,n.pca=PC_npcas,n.da=PC_ndims,scannf=FALSE) 
PC_test_a_score_all <- optim.a.score(PC_dapc_all)
PC_best_a_all = PC_test_a_score_all$best
PC_dapc_all <- dapc(PC_gpfile_all,PC_gpfile_all$pop,n.pca=PC_best_a_all,n.da=PC_ndims) 

PC_dapc_neu <- dapc(PC_gpfile_neu,PC_gpfile_neu$pop,n.pca=PC_npcas,n.da=PC_ndims,scannf=FALSE) 
PC_test_a_score_neu <- optim.a.score(PC_dapc_neu)
PC_best_a_neu = PC_test_a_score_neu$best
PC_dapc_neu <- dapc(PC_gpfile_neu,PC_gpfile_neu$pop,n.pca=PC_best_a_neu,n.da=PC_ndims) 

PC_dapc_out <- dapc(PC_gpfile_out,PC_gpfile_out$pop,n.pca=PC_npcas,n.da=PC_ndims,scannf=FALSE) 
PC_test_a_score_out <- optim.a.score(PC_dapc_out)
PC_best_a_out = PC_test_a_score_out$best
PC_dapc_out <- dapc(PC_gpfile_out,PC_gpfile_out$pop,n.pca=PC_best_a_out,n.da=PC_ndims) 

# make new df from loading coords & add pop and state cols
PC_dapc_all_df <- as.data.frame(PC_dapc_all$ind.coord)
PC_dapc_all_df$pop <- PC_popnames_gpfile 
PC_dapc_all_df$state <- PC_statenames_gpfile
PC_dapc_all_df$SS <- PC_reSS_gpfile
PC_dapc_all_df$NS <- PC_NSnames_gpfile

# make extra columns into factors
PC_dapc_all_df$pop <- factor(PC_dapc_all_df$pop,levels = PC_ordered_pops)
PC_dapc_all_df$state <- factor(PC_dapc_all_df$state,levels = PC_ordered_states)
PC_dapc_all_df$SS <- factor(PC_dapc_all_df$SS,levels = ordered_SS)
PC_dapc_all_df$NS <- factor(PC_dapc_all_df$NS,levels = ordered_NS)

# get string of percent pop var explained by first two axes
PC_perc_first_PC_all <- paste(substr(as.character((PC_dapc_all$eig[1] / sum(PC_dapc_all$eig))*100), 1, 5), "%", sep="")
PC_perc_second_PC_all <- paste(substr(as.character((PC_dapc_all$eig[2] / sum(PC_dapc_all$eig))*100), 1, 5), "%", sep="")

# make new df from loading coords & add pop and state cols
PC_dapc_neu_df <- as.data.frame(PC_dapc_neu$ind.coord)
PC_dapc_neu_df$pop <- PC_popnames_gpfile 
PC_dapc_neu_df$state <- PC_statenames_gpfile
PC_dapc_neu_df$SS <- PC_reSS_gpfile
PC_dapc_neu_df$NS <- PC_NSnames_gpfile

# make extra columns into  factors
PC_dapc_neu_df$pop <- factor(PC_dapc_neu_df$pop,levels = PC_ordered_pops)
PC_dapc_neu_df$state <- factor(PC_dapc_neu_df$state,levels = PC_ordered_states)
PC_dapc_neu_df$SS <- factor(PC_dapc_neu_df$SS,levels = ordered_SS)
PC_dapc_neu_df$NS <- factor(PC_dapc_neu_df$NS,levels = ordered_NS)

# get string of percent pop var explained by first two axes
PC_perc_first_PC_neu <- paste(substr(as.character((PC_dapc_neu$eig[1] / sum(PC_dapc_neu$eig))*100), 1, 5), "%", sep="")
PC_perc_second_PC_neu <- paste(substr(as.character((PC_dapc_neu$eig[2] / sum(PC_dapc_neu$eig))*100), 1, 5), "%", sep="")

# make new df from loading coords & add pop and state cols
PC_dapc_out_df <- as.data.frame(PC_dapc_out$ind.coord)
PC_dapc_out_df$pop <- PC_popnames_gpfile 
PC_dapc_out_df$state <- PC_statenames_gpfile
PC_dapc_out_df$SS <- PC_reSS_gpfile
PC_dapc_out_df$NS <- PC_NSnames_gpfile

# make extra columns into 
PC_dapc_out_df$pop <- factor(PC_dapc_out_df$pop,levels = PC_ordered_pops)
PC_dapc_out_df$state <- factor(PC_dapc_out_df$state,levels = PC_ordered_states)
PC_dapc_out_df$SS <- factor(PC_dapc_out_df$SS,levels = ordered_SS)
PC_dapc_out_df$NS <- factor(PC_dapc_out_df$NS,levels = ordered_NS)

# get string of percent pop var explained by first two axes
PC_perc_first_PC_out <- paste(substr(as.character((PC_dapc_out$eig[1] / sum(PC_dapc_out$eig))*100), 1, 5), "%", sep="")
PC_perc_second_PC_out <- paste(substr(as.character((PC_dapc_out$eig[2] / sum(PC_dapc_out$eig))*100), 1, 5), "%", sep="")

### DAPC by site

PC_DAPC_all_bypop <- ggplot(data = PC_dapc_all_df, aes(x= LD1, y= LD2)) + 
  geom_point(aes(colour = pop), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values=PC_site_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_perc_first_PC_all, y=PC_perc_second_PC_all) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  # # ggtitle("P. californicus") +
  # # theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="right") +
  theme(legend.title=element_blank())
PC_DAPC_all_bypop

PC_DAPC_neu_bypop <- ggplot(data = PC_dapc_neu_df, aes(x= LD1, y= LD2)) + 
  geom_point(aes(colour = pop), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values=PC_site_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_perc_first_PC_neu, y=PC_perc_second_PC_neu) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="right") +
  theme(legend.title=element_blank())
PC_DAPC_neu_bypop

PC_DAPC_out_bypop <- ggplot(data = PC_dapc_out_df, aes(x= LD1, y= LD2)) + 
  geom_point(aes(colour = pop), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values=PC_site_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_perc_first_PC_out, y=PC_perc_second_PC_out) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="right") +
  theme(legend.title=element_blank())
PC_DAPC_out_bypop

### DAPC by state

PC_DAPC_all_bystate <- ggplot(data = PC_dapc_all_df, aes(x= LD1, y= LD2)) + 
  geom_point(aes(colour = state), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values=PC_state_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_perc_first_PC_all, y=PC_perc_second_PC_all) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="right") +
  theme(legend.title=element_blank())
PC_DAPC_all_bystate

PC_DAPC_neu_bystate <- ggplot(data = PC_dapc_neu_df, aes(x= LD1, y= LD2)) + 
  geom_point(aes(colour = state), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values=PC_state_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_perc_first_PC_neu, y=PC_perc_second_PC_neu) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="right") +
  theme(legend.title=element_blank())
PC_DAPC_neu_bystate

PC_DAPC_out_bystate <- ggplot(data = PC_dapc_out_df, aes(x= LD1, y= LD2)) + 
  geom_point(aes(colour = state), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values=PC_state_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_perc_first_PC_out, y=PC_perc_second_PC_out) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="right") +
  theme(legend.title=element_blank())
PC_DAPC_out_bystate

### DAPC by SS

PC_DAPC_all_bySS <- ggplot(data = PC_dapc_all_df, aes(x= LD1, y= LD2)) + 
  geom_point(aes(colour = SS), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values=SS_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_perc_first_PC_all, y=PC_perc_second_PC_all) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="right") +
  theme(legend.title=element_blank())
PC_DAPC_all_bySS

PC_DAPC_neu_bySS <- ggplot(data = PC_dapc_neu_df, aes(x= LD1, y= LD2)) + 
  geom_point(aes(colour = SS), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values=SS_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_perc_first_PC_neu, y=PC_perc_second_PC_neu) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="right") +
  theme(legend.title=element_blank())
PC_DAPC_neu_bySS

PC_DAPC_out_bySS <- ggplot(data = PC_dapc_out_df, aes(x= LD1, y= LD2)) + 
  geom_point(aes(colour = SS), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values=SS_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_perc_first_PC_out, y=PC_perc_second_PC_out) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="right") +
  theme(legend.title=element_blank())
PC_DAPC_out_bySS

######################################################################
### Multi panel plots

### DAPC neutral v adaptive

PC_DAPC_neu_bypop <- ggplot(data = PC_dapc_neu_df, aes(x= LD1, y= LD2)) + 
  geom_point(aes(colour = pop), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values=PC_site_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_perc_first_PC_neu, y=PC_perc_second_PC_neu) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="right") +
  ggtitle("Putatively neutral loci") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title=element_blank())

PC_DAPC_neu_bypop

#PC_dapc_out_df$LD2 <- PC_dapc_out_df$LD2*-1
PC_DAPC_out_bypop <- ggplot(data = PC_dapc_out_df, aes(x= LD1, y= LD2)) + 
  geom_point(aes(colour = pop), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values=PC_site_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_perc_first_PC_out, y=PC_perc_second_PC_out) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="right") +
  ggtitle("Putatively adaptive loci") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title=element_blank())
PC_DAPC_out_bypop

ggarrange(PC_DAPC_neu_bypop, PC_DAPC_out_bypop, 
          ncol=2, nrow=1, 
          common.legend = TRUE, 
          legend="right")


### PCA neutral v adaptive
PC_gpfile_neu_df$Axis1 <- -1*PC_gpfile_neu_df$Axis1
PC_PCA_neu_bysite <- ggplot(data = PC_gpfile_neu_df, aes(x= Axis1, y= Axis2)) + 
  geom_point(aes(colour = Population), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values=PC_site_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_axis1_var_neu, y=PC_axis2_var_neu) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="right") +
  ggtitle("Putatively neutral") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title=element_blank())
PC_PCA_neu_bysite

PC_PCA_out_bysite <- ggplot(data = PC_gpfile_out_df, aes(x= Axis1, y= Axis2)) + 
  geom_point(aes(colour = Population), size = 3.0, alpha = 0.7) +
  scale_colour_manual(values=PC_site_colors) +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  labs(x=PC_axis1_var_out, y=PC_axis2_var_out) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="right") +
  ggtitle("Putatively neutral") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title=element_blank())
PC_PCA_out_bysite

ggarrange(PC_PCA_neu_bysite, PC_PCA_out_bysite, 
          ncol=2, nrow=1, 
          common.legend = TRUE, 
          legend="right")
