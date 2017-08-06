## 20170707 NL
## adapted from Katherine Silliman
## Code to subset one SNP per GBS locus from a VCF file. Chooses the SNP
## with the highest sample coverage. If there is a tie, chooses the 1st SNP in the loci. (may change to random)
## May be specific to VCF format output from ipyrad.

import sys
import linecache

# command line arguments: [1] input VCF, and [2] output VCF with one SNP per RAD tag
inputfile = sys.argv[1]
outputfile = sys.argv[2]

locidict = {}
lineNum = []
IN = open(inputfile, "r")
OUT = open(outputfile, "w")

# counts genotype line in file
n = 1

locus_list = []

for line in IN:
	if "#" not in line:
		linelist = line.strip().split()
		if "loc" in linelist[0]:
			loci = linelist[0]
		else:
			loci = int(linelist[0])
		#Column 8 is INFO column of VCF file
		NS = float(linelist[7].split(";")[0].split("=")[1])
		if loci not in locidict.keys():
			locidict[loci] = [NS,n]
			if loci not in locus_list:
				locus_list.append(loci)
		else:
			if locidict[loci][0] < NS:
				locidict[loci] = [NS,n]
	else:
		OUT.write(line)
	n += 1
IN.close()
print("Total SNPS: "+str(n)+"\nUnlinked SNPs: "+str(len(locidict.keys())))

for locus in locus_list:
	line = linecache.getline(inputfile, locidict[locus][1])
	OUT.write(line)
OUT.close()
