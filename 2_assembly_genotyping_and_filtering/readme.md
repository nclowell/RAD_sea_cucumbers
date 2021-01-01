#### Assembly, genotyping, filtering

This directory contains the notebooks and scripts associated with assembly, genotyping, and filtering.

First, I demultiplexed raw fastq files using the [*process_radtags*](https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php) (v2.4) from Stacks using the following parameters: -c, -r, -t 144, -q,  -s 0, --barcode_dist_1 2, -E phred33, and -e sbfI.

Second, I used the [*dDocent*](https://www.ddocent.com/) pipeline to assemble and genotype indivdiuals using default parameters. Configuration file [here](https://github.com/nclowell/SeaCukes/blob/master/2_assembly_genotyping_and_filtering/dDocent_config.txt).

Third, I filtered loci and individuals, following [this tutorial](https://www.ddocent.com/filtering/), using a combination of notebooks and scripts. To see how I filtered, start with [this notebook](https://github.com/nclowell/SeaCukes/blob/master/2_assembly_genotyping_and_filtering/filtering.ipynb).