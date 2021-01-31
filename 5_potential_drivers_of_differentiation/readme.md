## Drivers of population differentiation

I looked for evidence of isolation by distance, among all loci, putatively neutral loci, and putatively adaptive loci using this [script](https://github.com/nclowell/SeaCukes/blob/master/5_potential_drivers_of_differentiation/isolation_by_distance.R).

I also investigated biological functions associated with putatively adaptive SNPs to develop hypothese of selection driving population differentiation. I used the UniPort database to create a ``blastx`` database, and retrieved GO Slim terms for hits. Script [here](https://github.com/nclowell/SeaCukes/blob/master/5_potential_drivers_of_differentiation/gene_annotation_w_uniprot.R).

I used DAPC and PCA on putatively neutral versus putatively adaptive loci to look at potential differences in clustering patterns. That R script is located [here](https://github.com/nclowell/SeaCukes/blob/master/3_pop_structure_analyses/PCA_and_DAPC.R).

#### Some results

We found evidence for isolation by distance, with the majority of the correlation between geographic and genetic distance driven by putatively adaptive loci.

![img](https://github.com/nclowell/RAD_sea_cucumbers/blob/master/imgs/ibd.PNG?raw=true)



