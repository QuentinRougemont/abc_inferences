#!/usr/bin/env python
# script to generate input files for ABC
# ./script.py [data_genotypes] [species name 1] [species name 2] (optional [cutoff])
# cutoff = maximum percentage of missing data (per population) -> cutoff = 1 by default (no filtration)

print "--------------------------------------------------------------------------------"
print "script ABC Camille - modifie le 220816 par TL - thibault.leroy@pierroton.inra.fr"
print "--------------------------------------------------------------------------------"

#from colorama import *
import sys
print(sys.argv)

def treatLocus(locusSpecies1, locusSpecies2):
	locusSpecies1 = locusSpecies1.upper()
	locusSpecies2 = locusSpecies2.upper()
#	print lackingspecies1,lackingspecies2
	res = {}
	baseContent = {}
	res["nsp1"] = 0	# number of keeped individuals within species 1
	res["nsp2"] = 0 # number of keeped individuals within species 2
	res["seq1"] = ""	# sequence of species 1
	res["seq2"] = ""	# sequence of species 2
	res["statu"] = 0 # if == 0 -> the locus is rejected; if == 1 -> the locus is keeped for ABC analysis
	seqTot = locusSpecies1+locusSpecies2
	nsp1 = 0
	nsp2 = 0
	# count the number of possible caracters
	baseContent["A"] = seqTot.count("A")
	baseContent["T"] = seqTot.count("T")
	baseContent["C"] = seqTot.count("C")
	baseContent["G"] = seqTot.count("G")
	tmp = 0
	for i in baseContent.values():	# loop to count the number of absent nucleotide
		if i == 0:
			tmp += 1
	if tmp != 2:	# rejected locus if it doesn't contain a bi-allelic polymorphism
		res["nsp1"] = -9
		res["nsp2"] = -9
		res["seq1"] = "NA"
		res["seq2"] = "NA"
		res["statu"] = 0
		return(res)
	if tmp == 2:
		nN1 = locusSpecies1.count("-") + locusSpecies1.count("N")
		lengthseq1 = len(locusSpecies1)
		nN2 = locusSpecies2.count("-") + locusSpecies2.count("N")
		lengthseq2 = len(locusSpecies2)
		if nN1 != 0 or nN2 != 0: 
			locusSpecies1 = locusSpecies1.replace("-", "")
			locusSpecies1 = locusSpecies1.replace("N", "")
			locusSpecies2 = locusSpecies2.replace("-", "")
			locusSpecies2 = locusSpecies2.replace("N", "")
		if nN1 > cutoff * lengthseq1 or nN2 > cutoff * lengthseq2:
			print "discarded locus:", pouet, "because", nN1, " > ", cutoff, "*", lengthseq1, "or ", nN2, " > ", cutoff, "*", lengthseq2
			res["nsp1"] = -9
			res["nsp2"] = -9
			res["seq1"] = "NA"
			res["seq2"] = "NA"
			res["statu"] = 0
			return(res)
		else:
			res["nsp1"] = len(locusSpecies1)
			res["nsp2"] = len(locusSpecies2)
			res["statu"] = 1
			ancestralState = locusSpecies1[0]
			seq1 = "".join([ "0" if x is ancestralState else "1" for x in locusSpecies1 ])
			seq2 = "".join([ "0" if x is ancestralState else "1" for x in locusSpecies2 ])
			res["seq1"] = seq1
			res["seq2"] = seq2
			return(res)

inputFile = sys.argv[1]
species1 = sys.argv[2]
species2 = sys.argv[3]
cutoff = float(sys.argv[4]) if len(sys.argv) > 4 else 1 # cutoff = 1 by default (no filtration), if not set it by the 4th arguments
if cutoff == 1:
	print( 'CAUTION: your current parameters keep all loci, even those with too many missing data.')
#	print(Style.RESET_ALL)
else:
	print( 'current cutoff: '), cutoff
	print( 'All sites with more than [cutoff * number of haplotypes] per pop will be not be considered.')
#	print(Style.RESET_ALL)

# contains the indexes of columns for species 1 and species 2
colOfInterest = {}
colOfInterest[species1] = list()
colOfInterest[species2] = list()

file = open(inputFile, "r")

speciesNames = file.readline().strip().split("\t")
cnt = 0
for i in speciesNames:
	if i == species1:
		colOfInterest[species1].append(cnt)
	if i == species2:
		colOfInterest[species2].append(cnt)
	cnt += 1

# jump the useless second line
i = file.readline()

bpfileL1 = "# {0} {1} 1e+05".format(species1, species2)	# comment line
bpfileL2 = ""	# locus length (but useless for SNPs, so, equal to one)
bpfileL3 = ""	# nIndividuals for species 1
bpfileL4 = ""	# nIndividuals for species 2
bpfileL5 = ""	# contains theta for loci (-t argument), but equal to one for SNPs (-s argument)
bpfileL6 = ""	# contains rho for loci, but useless for SNPs

spinput = ""
msStyle = "./msnsam tbs 20 -s tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -ej tbs 1 2 -eN tbs tbs\n3579 27011 59243\n"

# loop over the loci
nlocus = 0
pouet = 0
for i in file:
	pouet += 1
	i = i.strip().split("\t")
	genotypesSpecies1 = ""
	genotypesSpecies2 = ""
	for j in colOfInterest[species1]:
		genotypesSpecies1 += i[j]
	for j in colOfInterest[species2]:
		genotypesSpecies2 += i[j]
	treated = treatLocus(genotypesSpecies1, genotypesSpecies2)
	if(treated["statu"] == 1):	# si on garde le locus
		nlocus += 1
		bpfileL2 += "{0}\t".format(1)
		bpfileL3 += "{0}\t".format(treated["nsp1"])
		bpfileL4 += "{0}\t".format(treated["nsp2"])
		bpfileL5 += "{0}\t".format(1)
		bpfileL6 += "{0}\t".format(1)
		msStyle += "\n//\t78\t1.58854\t1.58854\t111\t56\t22\t0.795745\t1.02097\t0.574558\t1.49299\t4.90703\t4.90703\t6.22109\nsegsites: 1\npositions:\t0.5\n"
		msStyle += "".join([ i+"\n" for i in treated["seq1"]+treated["seq2"] ])

		spinput += "{0}\n{1}\n{2}\n".format(treated["nsp1"], treated["nsp2"], 1)	# nsp1 // nsp2 // length of the locus (equal to one here)
file.close()


spinput = "\n{0}\n".format(nlocus) + spinput + "100000\nmyfifo\n"
bpfile = "{0}\n{1}\n{2}\n{3}\n".format(bpfileL1, bpfileL2, bpfileL3, bpfileL4)

outputFile = open("spinput.txt", "w")
outputFile.write(spinput)
outputFile.close()

outputFile = open("bpfile", "w")
outputFile.write(bpfile)
outputFile.close()

outputFile = open("locus.ms", "w")
outputFile.write(msStyle)
outputFile.close()

