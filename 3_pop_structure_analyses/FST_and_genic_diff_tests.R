#############################################################

# 20180628 NL
# Compute FST and run genic diff tests at pairwise and global scales

#############################################################

# import library
library("genepop")

# get arguments from command line
args <- commandArgs(TRUE)

# make objects out of arguments
gp_file = args[1]
Fst_outfile_pair = args[2]
Fst_outfile_glob = args[3]
Diff_outfile_pair = args[4]
Diff_outfile_glob = args[5]

# make a file object out of genepop file
infile <- system.file(gp_file, package="genepop")
locinfile <- gp_file
check <- file.copy(infile,locinfile,overwrite = TRUE)

# cal genepop's Fst and genic diff tests in genepop
Fst(locinfile, pairs = TRUE, outputFile = Fst_outfile_pair)
Fst(locinfile, pairs = FALSE, outputFile = Fst_outfile_glob)
test_diff(locinfile, genic = TRUE, pairs = TRUE, outputFile = Diff_outfile_pair)
test_diff(locinfile, genic = TRUE, pairs = TRUE, outputFile = Diff_outfile_glob)

