#### Population structure analyses

I calculated Fst and ran genic differentiation tests at the pairwise and global scales, using [this R script](https://github.com/nclowell/SeaCukes/blob/master/3_pop_structure_analyses/FST_and_genic_diff_tests.R).

I calculated expected heterozygosity, observed heterozygosity, and the proportion of polymorphic loci per site using [this R script](https://github.com/nclowell/SeaCukes/blob/master/3_pop_structure_analyses/get_He_Ho_propPolym_fromGP.R).

I plotted PCAs and DAPCs using [this R script](https://github.com/nclowell/SeaCukes/blob/master/3_pop_structure_analyses/PCA_and_DAPC.R).

I ran clustering analyses in ADMIXTURE following [this tutorial](https://speciationgenomics.github.io/ADMIXTURE/).

I used [this script]() to run AMOVAs using the following hierarchical groupings: by state/provine (arbitary mid-scale regions), inside or outside the Salish Sea, and North or South of the bifurcating North Pacific Current (for comparison with Xuereb et al. (2018))