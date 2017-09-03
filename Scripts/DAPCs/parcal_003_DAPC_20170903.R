# Let's first run a DAPC with all individuals and all loci


# Import these libraries
library(adegenet)
library(hierfstat)

# Set working directory
setwd("C:/Users/Natalie Lowell/SHARED_FOLDER/Learn_iPyrad/PARCAL_RUN1/parcal_004_outfiles/new_outfiles")


# Read in structure file, with a ".str" file extension
data_all_loci <-read.structure("parcal_004_biall_NOmaf_oneSNP_inames_pops.str")

names(data_all_loci)
data_all_loci$pop

BC_BellaBella <- rep("BC_BellaBella", 28)
AK_ChiniakBay <- rep("AK_ChiniakBay", 9) 
AK_AukeBay <- rep("AK_AukeBay", 19)
AK_YakutatBay <- rep("AK_YakutatBay", 21)

pop_groups <- as.factor(c(rep("BC_BellaBella", 28),rep("AK_ChiniakBay", 9), rep("AK_AukeBay",19), rep("AK_Yakutat",20)))                        
pop_labels <- c(BC_BellaBella, AK_ChiniakBay, AK_AukeBay, AK_YakutatBay)
pop_cols <- c("black","dodgerblue","tomato","deepskyblue")


# Data exploration to find # PCs. Go through all of their prompts:
dapc_all <- dapc(data_all_loci,data_all_loci$pop)

# [1] Choose the number PCs to retain (>=1):
#     You will see a plot of cumulative variance by number retained PCs
#     Retain N/3 Prinicpal components
#     E.g., if x axis goes until 70, retain 70/3 = 23.333 = 23
#     You can also get this number from the test_a_score plot, subtitle

# [2] Choose the number discriminant functions to retain (>=1):
#     You will see a fiery histogram with N number of bars
#     Retain the number of linear discriminants (I think that is == n.da)
#     

myclusters <- find.clusters(data_all_loci)
test_a_score <- optim.a.score(dapc_all) # this also gives you optimal number of PCs, subtitle of plot


dapc_all <- dapc(data_all_loci,data_all_loci$pop,n.pca=25,n.da=3) 

#2D plot
scatter(dapc_all,scree.da=FALSE,cellipse=0,leg=FALSE,label=c("BC_BellaBella","AK_ChiniakBay","AK_AukeBay","AK_Yakutat"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(15,15,15,15),solid=1)
legend(x = 2, y = -1,bty='n',legend=c("BC_BellaBella", "AK_ChiniakBay","AK_AukeBay","AK_Yakutat"),pch=c(15,15,15,15),col=pop_cols,cex=1.3)

