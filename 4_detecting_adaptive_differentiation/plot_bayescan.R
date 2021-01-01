#### 20191025 Natalie Lowell

### This code is used to explore the output of a bayescan analysis. 

##########################################################
# Load necassary libraries
library(ggplot2)
library(tidyverse)
library(genepopedit)
library(RColorBrewer)
library(cowplot)

plot_bayescan<-function(file,PO=-1,FDR=-1,w=-1,n=1000,size=1,pos=0.35)
{
  res=read.table(file)
  
  PO2FDR<-function(res,PO)
  {
    q=PO/(1+PO)
    significant=res[res$prob>=q,]
    FDR=sum(1-significant$prob)/nrow(significant)
    FDR
  }
  
  FDR2PO<-function(res,FDR)
  {
    q_vals=sort(res$prob)
    for (q in q_vals) {
      significant=res[res$prob>=q,]
      if(sum(1-significant$prob)/nrow(significant)<=FDR) {
        break
      }
    } 
    q/(1-q) 
  }
  
  
  BE<-function(res,PO,w)
  {
    q=PO/(1+PO)
    non_significant=res[res$prob<q,]
    FNDR=sum(non_significant$prob)/nrow(non_significant)
    
    significant=res[res$prob>=q,]
    FDR=sum(1-significant$prob)/nrow(significant)
    
    (1-w)*FDR+w*FNDR
  }
  
  if (FDR==-1 && PO>0 && w==-1) {
    FDR=PO2FDR(res,PO)  
  }
  else if (PO==-1 && FDR>0 && w==-1) {
    PO=FDR2PO(res,FDR)
    FDR=PO2FDR(res,PO)
  }
  else if (PO==-1 && FDR==-1 && w>=0) {
    logPO=seq(0,4,length=n)
    errors=0
    for (i in 1:n)
    {
      errors[i]=BE(res,10^(logPO[i]),w)
    }
    PO=10^(logPO[max(which(errors==min(errors)))])
    FDR=PO2FDR(res,PO)
  }
  else {  
    stop("One of PO, FDR or w is required as an argument of the function")
  }
  if (PO==Inf)
    p=1
  else
    p=PO/(1+PO)
  non_significant=res[res$prob<p,]
  FNDR=sum(non_significant$prob)/nrow(non_significant)
  
  # when p=1 replaces by p=0.9999 to plot with logit scale
  if(nrow(res[res$prob>=0.9999,])>=1)
    res[res$prob>=0.9999,]$prob=0.9999
  if(nrow(res[res$prob<=0.001,])>=1)
    res[res$prob<=0.001,]$prob=0.001
  
  outliers=as.integer(row.names(res[res$prob>=p,]))
  
  ok_outliers=TRUE
  if (PO==Inf | sum(res$prob>=p)==0)
    ok_outliers=FALSE;
  
  if (w>=0) {
    #par(mfrow=c(2,1))
    par(mar = c(5, 4, 4, 4) + 0.3)
  }
  # plot
  plot(log10(res$prob/(1-res$prob)),res$fst,xlab="log(PO)",ylab="Fst",pch=19,cex=size)
  # add names of loci over p and vertical line
  if (ok_outliers) {
    text(log10(res[res$prob>=p,]$prob/(1-res[res$prob>=p,]$prob))+pos*(round(runif(nrow(res[res$prob>=p,]),1,2))*2-3),res[res$prob>=p,]$fst,row.names(res[res$prob>=p,]),cex=size)
    lines(c(log10(PO),log10(PO)),c(-1,1),lwd=2)
  }
  
  if (w>=0) {
    par(new = TRUE)
    plot(logPO,errors,xlim=range(log10(res$prob/(1-res$prob))), type = "l", axes = FALSE, lty=3, bty = "n", xlab = "", ylab = "")
    #  plot(logPO,errors,type="l",xlab="log(PO)",ylab="Posterior Odds")
    axis(side=4, at = pretty(range(errors)))
    mtext("Posterior Odds", side=4, line=3)
  }
  #par(mfrow=c(1,1))
  return(list("PO"=PO,"FDR"=FDR,"FNDR"=FNDR,"p"=p,"outliers"=outliers,"nb_outliers"=length(outliers)))
}

###########################################################

setwd("E:/dDocent_for_mox/parcal_wd/parcal_mox001/outliers")

# Specify the file names for the input files
bayescan_file <- "parcal_mox001_bscan_out_Fst.txt" 
genepop_file <- "primSNPs_noINDL_parcal_mox001_md70_maf05_minQ20_mDP10_inames_noreps_HWE_oneSNPhiMAF_gorder.txt"

##########################################################
# Specify the file names for the output files

bayescan_output_df_all <- "bayescan_parcalmox001_all_loci.txt"
bayescan_output_df_outliers <- "bayescan_significants_parcalmox001_all.txt"

################################################################
################################################################
# Part 1 : Read in the bayescan results and plot some summary stats

mydata <- read.table(bayescan_file, header = TRUE)
head(mydata)
length(mydata$bayescan_index)

# Plot some summary statistics from the bayescan results

# plot the distribution of fst

plot_fst <- ggplot() +
  geom_histogram(data =  mydata, aes(fst), bins = 100)+
  xlab(expression(italic(F[ST])))+
  ylab("Number of loci") +
  theme_classic()+
  geom_hline(yintercept=0, colour="white", size=0.25) # add a fake white line to remove the obnoxious "tracer line" at bottom of plot

plot_fst

plot_fsttail <- ggplot() +
  geom_histogram(data =  mydata, aes(fst), bins = 100)+
  xlab(expression(italic(F[ST])))+
  ylab("Number of loci") +
  coord_cartesian(xlim = c(0, 0.26), ylim = c(0, 100))+
  theme_classic()

plot_fsttail

# plot the distribution of q values

plot_qval <- ggplot() +
  geom_histogram(data =  mydata, aes(qval), bins = 100)+
  xlab("q-value") +
  ylab("Number of loci")+
  theme_classic()

plot_qval

# plot the fst vs qvalue

ggplot()+
  geom_point(data = mydata, aes(x = qval, y = fst), shape = 1, alpha = 0.4)+
  ylab(expression(italic(F[ST]))) +
  xlab("q-value")+
  theme_classic()


# plot the fst vs prob

ggplot()+
  geom_point(data = mydata, aes(x = prob, y = fst), shape = 1, alpha = 0.4)+
  ylab(expression(italic(F[ST]))) +
  xlab("Posterior probability")

# plot the fst vs log10_prob

plot_log10 <- ggplot()+
  geom_point(data = mydata, aes(x = log10.PO., y = fst), shape = 1, alpha = 0.4)+
  ylab(expression(italic(F[ST])))+
  xlab("log10 (posterior probability)")+
  theme_classic()

plot_log10

# arrange the three plots in a single row
multiplot <- plot_grid( plot_fst,
                        plot_fsttail ,
                        plot_qval, 
                        plot_log10 ,
                        align = 'vh',
                        labels = c("A", "B", "C", "D"),
                        hjust = 0,
                        nrow = 2)
multiplot

### Get outlier locus names

mydata <- read.table(bayescan_file, header = TRUE)

# read in locus names
locus_names <- read.csv("../oneSNPhiMAF_2075L_names.txt", header=FALSE)

# add names column
data_inames <- cbind(locus_names, mydata)

# pick a threshold!! Here... 750
threshold = 750

# subset df for rows where log 10 posterior prob is over threshold
outliers_df <- filter(data_inames, log10.PO. > threshold)
outliers_df

# write to file
write.table(outliers_df, file = "bayescan_outliers.txt", sep="\t")

p_b_results <- plot_bayescan(bayescan_file,FDR=0.05)
bs_out_indeces <- p_b_results$outliers
length(bs_out_indeces)

mydata2 <- mydata %>%
  mutate(BS_out_FDR=FALSE) %>%
  rowwise() %>%
  mutate(log10_q = log10(qval))
for(idx in bs_out_indeces){
  mydata2[idx,"BS_out_FDR"] = TRUE
}
data2_inames <- cbind(locus_names, mydata2)

sig_data2_inam <- data2_inames %>%
  filter(BS_out_FDR == TRUE)

write_csv(sig_data2_inam, "Bayescan_50outliers_FDR_20201013.csv")

plot_log10q_fst <- ggplot()+
  geom_point(data = mydata2, aes(x = log10_q, 
                                 y = fst, color=BS_out_FDR), 
             alpha = 0.4)+
  ylab(expression(italic(F[ST]))) +
  scale_color_manual(values=c("black", "red")) +
  xlab("log10 q-value") +
  
  theme_classic() +
  #scale_x_reverse() +
  xlim(0,-6) +
  theme(legend.position="none") #+

plot_log10q_fst





