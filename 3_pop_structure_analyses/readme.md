#### Population structure analyses

I calculated Fst and ran genic differentiation tests at the pairwise and global scales, using [this R script](https://github.com/nclowell/SeaCukes/blob/master/3_pop_structure_analyses/FST_and_genic_diff_tests.R).

I calculated expected heterozygosity, observed heterozygosity, and the proportion of polymorphic loci per site using [this R script](https://github.com/nclowell/SeaCukes/blob/master/3_pop_structure_analyses/get_He_Ho_propPolym_fromGP.R).

I plotted PCAs and DAPCs using [this R script](https://github.com/nclowell/SeaCukes/blob/master/3_pop_structure_analyses/PCA_and_DAPC.R).

I ran clustering analyses in ADMIXTURE following [this tutorial](https://speciationgenomics.github.io/ADMIXTURE/).

I used [this script](https://github.com/nclowell/RAD_sea_cucumbers/blob/master/3_pop_structure_analyses/amova.R) to run AMOVAs using the following hierarchical groupings: by state/provine (arbitary mid-scale regions), inside or outside the Salish Sea, and North or South of the bifurcating North Pacific Current (for comparison with Xuereb et al. (2018))

#### Some population structure results

Overall, detectable population structure almost every way we sliced it. Global Fst = 0.0068 and pairwise Fsts are presented in the table below. Asterisks denote signficant genic differentiation tests.

![fst](https://github.com/nclowell/RAD_sea_cucumbers/blob/master/imgs/pairwise_fst_table.PNG?raw=true)