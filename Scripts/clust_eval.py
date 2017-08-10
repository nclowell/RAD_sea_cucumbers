###### Clustering similarity evaluation + plots #######
# PURPOSE: make plots to help user evaluate effect of clustering
# similarity parameter value on genotypes in ipyrad
# INPUTS: managed by argparse, include:
# - path to ipyrad directory
# - text file with all assembly names, each on own line
# - text file with all clustering param values, each on
#   own line and matching order of assembly names text file
# OUTPUTS: plots comparing many output metrics to help evaluate
# the effect of clustering similarity
######################################################

# import modules
from __future__ import division
import os
import argparse
import matplotlib.pyplot as plt
import numpy as np

# manage command line arguments
parser = argparse.ArgumentParser(description="Read in step 7 stats files, plot contents to help user evaluate effect of clustering similarity")
parser.add_argument("-d", "--dir", help="path to ipyrad directory", type=str, required = True)
parser.add_argument("-n", "--names", help="text file of assembly names, each on own line, stored in ipyrad dir", type=str, required = True)
parser.add_argument("-c", "--clusts", help="text file of clustering similarity parameter values, each on own line, and matching order of assembly names text file, stored in ipyrad dir", type=str, required = True)
args = parser.parse_args()

# change working directory to ipyrad directory
os.chdir(args.dir)

# get assembly names out of input file
assemblies_file = open(args.names,"r")
assemblies = [] #initiate list to store assembly names
for line in assemblies_file:
    name = line.strip()
    if name != "":
        assemblies.append(name)

# get clustering similarity parameter values out of input file
clusts_file = open(args.clusts,"r")
clusts = []
for line in clusts_file:
    clust = line.strip()
    if clust != "":
        clusts.append(float(clust))

# initiate lists to store values of different metrics
pi_0 = []
var_0 = []
total_filters_MA = []
total_filters_RD = []
tot_filt_loci = []
tot_prefilt_loci = []
ret_loci_after_minsample = []
applied_order_MA = []

# iterate through assembles, extract information from stats files, and store in lists
for assembly in assemblies:
    stats_file = open(assembly + "_outfiles/" + assembly + "_stats.txt", "r")
    stats_file_lines = stats_file.readlines()
    total_filters_MA.append(float(stats_file_lines[12].strip().split()[1]))
    total_filters_RD.append(float(stats_file_lines[7].strip().split()[1]))
    tot_filt_loci.append(float(stats_file_lines[13].strip().split()[1]))
    tot_prefilt_loci.append(float(stats_file_lines[6].strip().split()[1]))
    ret_loci_after_minsample.append(float(stats_file_lines[11].strip().split()[3]))
    applied_order_MA.append(float(stats_file_lines[12].strip().split()[2]))
    for line in stats_file_lines:
        if line[0:6] == "## pis":
            index = stats_file_lines.index(line)
            pi_0.append(float(stats_file_lines[index+4].strip().split()[3]))
            var_0.append(float(stats_file_lines[index+4].strip().split()[1]))

### MAKE PLOTS

# make dir for plots
if not os.path.exists("clust_eval_plots"):
    os.makedirs("clust_eval_plots")

# proportion homozygous loci with singletons
plt.plot(clusts,[x/y for x, y in zip(var_0, tot_filt_loci)])
plt.title("Proportion of homozygous loci including singletons")
plt.xlabel("Clustering similarity %")
plt.ylabel("var = 0 / total filtered loci")
plt.savefig("clust_eval_plots/prop_hom_wsingletons.png")
plt.close()

# proportion homozygous loci without singletons
plt.plot(clusts,[x/y for x, y in zip(pi_0, ret_loci_after_minsample)])
plt.title("Proportion of loci homozygous without singletons")
plt.xlabel("Clustering similarity %")
plt.ylabel("pi = 0 / retained loci after minimum sample")
plt.savefig("clust_eval_plots/prop_hom_wosingletons.png")
plt.close()

# proportion heterozygous loci
props_het = []
for i in [(y-x)/y for x, y in zip(pi_0, tot_filt_loci)]:
    props_het.append(1-i)
plt.plot(clusts,props_het)
plt.title("Proportion of heterozygous loci")
plt.xlabel("Clustering similarity %")
plt.ylabel("(total filtered loci - pi_0)/total filtered loci")
plt.savefig("clust_eval_plots/prop_het.png")
plt.close()

# proportion filtered out for polyploidy (max alleles = 2)
plt.plot(clusts,[x/y for x, y in zip(applied_order_MA, ret_loci_after_minsample)])
plt.title("Proportion loci filtered out for polyploidy")
plt.xlabel("Clustering similarity %")
plt.ylabel("Applied_order filtered max alleles / retained_loci after min sample")
plt.savefig("clust_eval_plots/prop_lost_polyploids.png")
plt.close()

# proportion loci lost to remove duplicates (potentially over clustering)
plt.plot(clusts,[x/y for x, y in zip(total_filters_RD, tot_prefilt_loci)])
plt.title("Proportion of loci lost to remove_duplicates")
plt.xlabel("Clustering similarity %")
plt.ylabel("total_filters remove duplicates / total prefiltered loci")
plt.savefig("clust_eval_plots/prop_lost_remdups.png")
plt.close()

# 20170810 Natalie Lowell
