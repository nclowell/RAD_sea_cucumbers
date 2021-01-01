#### Detecting putatively adaptive differentiation

We used several methods for detected putatively adaptive differentiation.

The first is Fst outlier detection methods, for which we used [*BayeScan*](http://cmpg.unibe.ch/software/BayeScan/) and [*OutFLANK*](http://rstudio-pubs-static.s3.amazonaws.com/305384_9aee1c1046394fb9bd8e449453d72847.html).

The second is gene environment association methods, for which we used [*bayenv2*](https://gcbias.org/bayenv/) and redundancy analysis. For both, I gathered environmental predictor data from the Bio-Oracle and Bio-Oracle2 databases using [this script](https://github.com/nclowell/SeaCukes/blob/master/4_detecting_adaptive_differentiation/access_biooracle_for_env_predictors.R). Environmental data was not available for one site (Eld Inlet, WA) so it was excluded from further gene-environment association analyses. 

Then, for bayenv2, I ran it using scaled environmental data, using default parameters (-k 100000) and the supplied script for running *bayenv2* per locus. Results were parsed using this script.