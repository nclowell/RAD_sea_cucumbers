#################### CONVERT GENEPOP TO OUTFLANK INPUTS #######################

## Purpose: This script will produce the required input for OutFlank's function 
## MakeDiploidFSTMat(), which creates the required input file format for 
## OutFLANK. 
##
## Inputs:
## - Genepop file with a header line, locus names one per line
## - population map file, as in ipyrad or Stacks
##
## Outputs:
## - SNP matrix file
## - locus names file, one per line matching order in genepop
## - pop names file, one per line matching order in genepop
##
## Assumptions:
## - ind nor locus names begin with "pop" or "Pop"
##

###############################################################################


import argparse
import numpy as np 

parser = argparse.ArgumentParser(description="produce SNPmat file, and files containing loci / population lists for OutFLANK outlier analysis.")

parser.add_argument("-i", "--input", help="genepop file that you want to run through OutFLANK")
parser.add_argument("-p", "--popmap", help="population map from stacks (each line has sample - tab - population")
parser.add_argument("-s", "--snpmat", help="snp matrix file for OutFLANK")
parser.add_argument("-ol", "--outLocusNames", help="text file with the name of each locus on each line, to be read into R")
parser.add_argument("-op", "--outPopNames", help="text file with the name of each sample's population on each line, to be read into R")

args = parser.parse_args()

infile = open(args.input, "r")
popmap = open(args.popmap, "r")
snpmat = open(args.snpmat, "w")
locusfile = open(args.outLocusNames, "w")
popfile = open(args.outPopNames, "w")


########### Make SNPmat file and get locusfile output #########

# make locus name input file for OutFLANK
loci_str = ""
lines = infile.readlines()
if lines[0][0] == "#":
    print("Header of genepop file: " + "\n" + lines[0])
else:
    loci_str += '"' + lines[0].strip() + '"' + "\n"

check = 0	
for line in lines[1:]:
    if line.startswith("Pop") or line.startswith("pop"):
        check = 1
    if check == 0:
        loci_str += '"' + line.strip() + '"' + "\n"
locusfile.write(loci_str.strip())
locusfile.close()

# get array of genotypes, transpose
gen_list_for_array = []
for line in lines:
    temp_gen_list = []
    if line[0] != "#" and len(line.strip().split()) > 2: # find genotype lines
        broken_line = line.split(",")
        for genotype in broken_line[1].strip().split():
            temp_gen_list.append(genotype)
        gen_list_for_array.append(temp_gen_list)
gen_array = np.array(gen_list_for_array)
gen_tarray = np.transpose(gen_array)	
    
snp_matrix_list_for_array = []
for row in gen_tarray:
    
    # store allele counts
    allele_dict = {}
    for genotype in row:
        allele1 = genotype[:int(len(genotype)/2)]
        allele2 = genotype[int(len(genotype)/2):]
        if allele1 not in allele_dict:
            allele_dict[allele1] = 1
        else:
            allele_dict[allele1] += 1
        if allele2 not in allele_dict:
            allele_dict[allele2] = 1
        else:
            allele_dict[allele2] += 1
    
    # find major allele
    max_count = 0
    for allele in allele_dict.keys():
        if allele_dict[allele] > max_count:
            major_allele = allele
            max_count = allele_dict[allele]
    
    snp_matrix_row_list = []
    for genotype in row:
        if float(genotype) == 0:
            snp_matrix_row_list.append(9)
        elif genotype == major_allele+major_allele:
            # homozygote for major allele
            snp_matrix_row_list.append(2)
        elif genotype[:int(len(genotype)/2)] == genotype[int(len(genotype)/2):]:
            # homozygote for minor allele
            snp_matrix_row_list.append(0)
        else:
            # heterozygote
            snp_matrix_row_list.append(1)
    if len(snp_matrix_row_list) > 0 :
        snp_matrix_list_for_array.append(snp_matrix_row_list)
   

snp_matrix_array = np.array(snp_matrix_list_for_array)
snp_matrix_tarray = np.transpose(snp_matrix_array)

for row in snp_matrix_tarray:
    row_str = ""
    for geno in row:
        row_str += str(geno) + " "
    row_str = row_str[:-1]
    snpmat.write(row_str + "\n")
snpmat.close()	

####### Make Popfile output #########

pop_dict = {}
for line in popmap:
	pop_dict[line.strip().split()[0]] = line.strip().split()[1]
popmap.close()
infile = open(args.input, "r")
infile.readline() # header, may have len(line.strip().split()) > 1
pops_str = ""
for line in infile:
	linelist = line.strip().split()
	if len(linelist) > 1:
		sample = linelist[0].strip(",")
		pops_str += '"' + pop_dict[sample] + '"\n'
infile.close()
popfile.write(pops_str.strip())
popfile.close()
