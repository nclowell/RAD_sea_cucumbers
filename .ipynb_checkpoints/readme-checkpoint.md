## Population structure and adaptive differentiation in the Giant California sea cucumber

![img](https://github.com/nclowell/SeaCukes/blob/master/Imgs_for_Notebooks/cukepic.jpg?raw=true)
<br>Photo credit: Pacific Shellfish Institute


#### Background
A research team led by Dan Cheney of the Pacific Shellfish Institute is developing methods to co-culture the giant California sea cucumber, *Parastichopus californicus*, underneath currently farmed species like salmonids, sable fish, and oysters. One part of this study is to investigate whether spatially distinct populations of the species exist, as to inform out-planting strategies. To answer this question, I am investigating whether there is population structure in *Parastichopus californicus* using RAD sequencing. 

This research was funded by NOAA Saltonstall-Kennedy program grant number 2004246. A one-page project description can be found [here](http://www.pacshell.org/pdf/SK%20Sea%20Cuc.pdf).

#### Project Goals

**[1]** Quantify and characterize patterns in population structure
<br>**[2]** Identify putatively adaptive loci and evaluate whether patterns of population structure differ between putatively neutral and adaptive data sets

#### Overview of methods

I extracted, quantified and verified quality of DNA from 9 sites spanning from Alaska to Oregon, and then prepared single-digest RAD libraries using *sbfI* follwing [Etter et al 2011](https://link.springer.com/protocol/10.1007/978-1-61779-228-1_9). DNA was sequenced to 150bp at BGI. You can find the library prep protocols in [1_library_prep](https://github.com/nclowell/SeaCukes/tree/master/1_library_prep).

I clustered loci and genotyped individuals using the [*dDocent*](https://www.ddocent.com/) pipeline, using the reference genome of a closely related species, *Parastichopus parvimensis*. Genotype data were filtered using [*vcftools*](http://vcftools.sourceforge.net/) and custom python scripts, such that retained loci passed the following filters:

- minimum minor allele count of 5 reads
- minimum quality score of 20
- minimum genotype depth of 10 reads
- minimum minor allele frequency of 5%
- maximum missing data per locus of 30% across sites
- one SNP per RAD tag, retaining that with the highest minor allele frequency
- in Hardy Weinberg Equilibrium

Individuals were retained if they were successfully genotyped at at least 70% of loci. You can find documentation for my assembly, genotyping, and filtering [here](https://github.com/nclowell/SeaCukes/tree/master/2_assembly_genotyping_and_filtering).

We used [*genepop*](https://cran.r-project.org/web/packages/genepop/index.html) to calculate FST and run exact G tests at global and pairwise scales. We used *genepop* to calculate expected and observed heterozygosity as well, and used custom python scripts to calculate the proportion of polymorphic loci per site.

We used principal component analysis and dicriminant analysis of principal components to investigate patterns of population differentiation, using the R package [*adegenet*](https://cran.r-project.org/web/packages/adegenet/index.html). We used [*ADMIXTURE*](https://gaworkshop.readthedocs.io/en/latest/contents/07_admixture/admixture.html) to conduct a clustering analysis and the R package [*poppr*](https://cran.r-project.org/web/packages/poppr/index.html) to summarize the partitioning of variance among different hierarchical groupings in AMOVAs. We tested for isolation by distance using a Mantel test in R.

Loci were classified as putatively adaptive if identified using at least one of two approaches: FST outlier detection and gene-environment association. We used [*Bayescan*](http://cmpg.unibe.ch/software/BayeScan/) and the R package [*OutFLANK*](http://rstudio-pubs-static.s3.amazonaws.com/305384_9aee1c1046394fb9bd8e449453d72847.html) to detect FST outliers. We used [*bayenv2*](https://gcbias.org/bayenv/) and redundancy analysis to test for gene-environment association and identify SNPs putatively involved in local adaptation. We used [*blastx*](https://blast.ncbi.nlm.nih.gov/Blast.cgi?LINK_LOC=blasthome&PAGE_TYPE=BlastSearch&PROGRAM=blastx) and the [*UniProt Knowledge Base*](https://www.uniprot.org/help/uniprotkb) to identify potentially associated biological processes for putatively adaptive loci. 

Once putatively adaptive SNPs were identified, they were used to distinguish a putatively neutral data set (by excluding these SNPs) and a putatively adaptive data set (the remaining SNPs). These were used to generate hypotheses about the roles of neutral versus adaptive forces in driving population differentiation.

#### This repository

This repository contains the protocols, scripts, and notebooks that I used in this project. Here are the quick links:

1. [1_library_prep](https://github.com/nclowell/SeaCukes/tree/master/1_library_prep)
2. [2_assembly_genotyping_and_filtering](https://github.com/nclowell/SeaCukes/tree/master/2_assembly_genotyping_and_filtering)


