#!/usr/bin/env python
# -*- coding: utf-8 -*-
from numpy.random import uniform
from numpy.random import beta
from numpy.random import binomial
from random import shuffle

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 

import sys
n1, n2, nA, tau, M1, M2, shape1, shape2 = [], [], [], [], [], [], [], []
help="\n\t\033[32mExternal required library: numpy \033[1;m(sudo apt-get install python-numpy)\n\t\
priorgen.py generates prior distributions for multiple multilocus simulations under 14 different models of speciation.\
The output can be used from the stdout by ms (Hudson 2002), msnsam (Ross-Ibarra 2008) and msms (Ewing and Hermisson 2010) using the 'tbs' feature.\n\tIt requires one input file containing six lines: \n\
\t\tL1=description line, non read by priorgen.py\n\
\t\tL2=a vector with the number of SNPs to simulate \033[1;35m(For Oak's data, this line is useless, so, read but neglected)\033[1;m\n\
\t\tL3=a vector with the number of sampled individuals (nspA) for each locus, for the first population\n\
\t\tL4=a vector with the number of sampled individuals (nspB) for each locus, for the second population\n\
\tValues print in the stdout are used by ms-like coalescent simulators, values written in a file are the multilocus parameters useful for an ABC analysis\n\n\
\t\tparameters: name of the output file name. Ex \033[1;35mparameters=listOfParameters.txt\033[1;m\n\
\t\tn1: prior for N1 (the effective population size of the first population). Ex \033[1;35mn1=0 n1=10\033[1;m\n\
\t\tn2: prior for N2 (the effective population size of the second population). Ex\033[1;35m n2=0 n2=10\033[1;m\n\
\t\tnA: prior for NA (the effective population size of the ancetral population). Ex\033[1;35m nA=0 nA=10\033[1;m\n\
\t\ttau: prior for Tsplit (the time of speciation). Ex\033[1;35m tau=0 tau=3\033[1;m\n\
\t\tM1 (M2): prior for migration rate 4.N1.m1 (4.N2.m2) into the first (second) population. Ex\033[1;35m M1=0 M1=4 M2=0 M2=4\033[1;m\n\
\t\tshape1: prior for the first shape parameter of the Beta distribution. Ex\033[1;35m shape1=0 shape1=10\033[1;m\n\
\t\tshape2: prior for the second shape parameter of the Beta distribution. Ex\033[1;35m shape2=0 shape2=50\033[1;m\n\
\t\tmodel: \033[1;35m=SI\033[1;m (Strict Isolation), \033[1;35m=IM\033[1;m (Isolation with Migration), \033[1;35m=AM\033[1;m (Ancient Migration) or \033[1;35m=SC\033[1;m (Secondary Contact)\n\
\t\tNvariation: \033[1;35m=homo\033[1;m (shared values of Ne throughout genome for N1, N2 and Nanc) or \033[1;35m=hetero\033[1;m (variation of Ne throughout genome for N1, N2 and Nanc)\n\
\t\tMvariation: \033[1;35m=homo\033[1;m (shared values of M throughout genome) or \033[1;35m=hetero\033[1;m (variation of M throughout genome)\n\
\t\tnreps: number of multilocus simulations. Ex \033[1;35mnreps=1000\033[1;m\n\
\t\tsymMig: \033[1;35m=sym\033[1;m (M1=M2) or \033[1;35m=asym\033[1;m (M1 and M2 are independently chosen)\n\n\
\tEx:\n\
\t\033[1;32m./priorgen.py bpfile=bpfile_test.txt n1=0 n1=1 n2=1 n2=2 nA=2 nA=3 tau=3 tau=4 M1=4 M1=5 M2=5 M2=6 shape1=0 shape1=10 shape2=0 shape2=100 model=IM nreps=2 Nvariation=homo Mvariation=hetero symMig=asym parameters=output.txt\033[1;m\n\n\
\tMore details about the models in \033[33mhttp://onlinelibrary.wiley.com/doi/10.1111/jeb.12425/abstract\033[1;m and \033[33mhttp://mbe.oxfordjournals.org/content/30/7/1574\033[1;m\n\n\
\tmsnsam tbs 20000 -s 1 -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -ej tbs 2 1 -eN tbs tbs\t#for 'SI'\n\
\tmsnsam tbs 20000 -s 1 -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -ej tbs 2 1 -eN tbs tbs\t#for 'IM'\n\
\tmsnsam tbs 20000 -s 1 -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs\t#for 'AM'\n\
\tmsnsam tbs 20000 -s 1 -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs\t#for 'SC'\n\n\
\tcamille.roux.1@unil.ch\n\t26/04/2015\n"
for i in sys.argv:
	if("help" in i):
		print(help)
		sys.exit(0)

if len(sys.argv)<=1:
	print(help)
	sys.exit(0)

for i in sys.argv:
	if "=" in i:
		i=i.split("=")
		if(i[0]=="bpfile"):
			bpfile=i[1]
		if(i[0]=="parameters"):
			outputParameters=i[1]
		if(i[0]=="n1"):
			n1.append(float(i[1]))
		if(i[0]=="n2"):
			n2.append(float(i[1]))
		if(i[0]=="nA"):
			nA.append(float(i[1]))
		if(i[0]=="tau"):
			tau.append(float(i[1]))
		if(i[0]=="M1"):
			M1.append(float(i[1]))
		if(i[0]=="M2"):
			M2.append(float(i[1]))
		if(i[0]=="shape1"):
			shape1.append(float(i[1]))
		if(i[0]=="shape2"):
			shape2.append(float(i[1]))
		if(i[0]=="model"):
			model=i[1]
		if(i[0]=="nreps"):
			nreps=int(i[1])
		if(i[0]=="Nvariation"):
			Nvariation=i[1]
		if(i[0]=="Mvariation"):
			Mvariation=i[1]
		if(i[0]=="symMig"):
			sym=i[1]

def binomBeta(nlocus, shape1, shape2, scalar):
	neutre=[0]
	hetero=[1]
	nNeutre=int(uniform(0, nlocus))
	nHetero=nlocus-nNeutre
	status=nNeutre*neutre+nHetero*hetero
	shuffle(status)
	values=[]
	for i in status:
		if i==0:
			values.append(scalar)
		if i==1:
			values.append(scalar*beta(shape1, shape2))
	res={}
	res["values"]=values
	res["nNeutre"]=nNeutre
	return(res)

infile=open(bpfile, "r")
tmp=infile.readline()	#skip the header
#L=infile.readline().strip().replace(" ", "\t")
#L=L.split("\t")
S = infile.readline().strip().replace(" ", "\t")
S = S.split("\t")
nspA = infile.readline().strip().replace(" ", "\t")
nspA = nspA.split("\t")
nspB = infile.readline().strip().replace(" ", "\t")
nspB = nspB.split("\t")
#theta=infile.readline().strip().replace(" ", "\t")
#theta=theta.split("\t")
#rho=infile.readline().strip().replace(" ", "\t")
#rho=rho.split("\t")

nlocus=int(len(S))

if(model=="SI"):
	if(Nvariation=="homo"):
		res="N1\tN2\tNa\tTsplit\n"
	if(Nvariation=="hetero"):
		res="N1\tN2\tNa\tTsplit\tshape1Ne\tshape2Ne\tpropNtrlNe1\tpropNtrlNe2\n"
if(model=="IM"):
	if(Nvariation=="homo"):
		if(Mvariation=="homo"):
			res="N1\tN2\tNa\tTsplit\tM1\tM2\n"
		if(Mvariation=="hetero"):
			res="N1\tN2\tNa\tTsplit\tM1\tM2\tshape1M1\tshape2M1\tshape1M2\tshape2M2\tpropNtrlM1\tpropNtrlM2\n"
	if(Nvariation=="hetero"):
		if(Mvariation=="homo"):
			res="N1\tN2\tNa\tTsplit\tshape1Ne\tshape2Ne\tpropNtrlNe1\tpropNtrlNe2\tM1\tM2\n"
		if(Mvariation=="hetero"):
			res="N1\tN2\tNa\tTsplit\tshape1Ne\tshape2Ne\tpropNtrlNe1\tpropNtrlNe2\tM1\tM2\tshape1M1\tshape2M1\tshape1M2\tshape2M2\tpropNtrlM1\tpropNtrlM2\n"
if(model=="AM"):
	if(Nvariation=="homo"):
		if(Mvariation=="homo"):
			res="N1\tN2\tNa\tTsplit\tTam\tM1\tM2\n"
		if(Mvariation=="hetero"):
			res="N1\tN2\tNa\tTsplit\tTam\tM1\tM2\tshape1M1\tshape2M1\tshape1M2\tshape2M2\tpropNtrlM1\tpropNtrlM2\n"
	if(Nvariation=="hetero"):
		if(Mvariation=="homo"):
			res="N1\tN2\tNa\tTsplit\tTam\tshape1Ne\tshape2Ne\tpropNtrlNe1\tpropNtrlNe2\tM1\tM2\n"
		if(Mvariation=="hetero"):
			res="N1\tN2\tNa\tTsplit\tTam\tshape1Ne\tshape2Ne\tpropNtrlNe1\tpropNtrlNe2\tM1\tM2\tshape1M1\tshape2M1\tshape1M2\tshape2M2\tpropNtrlM1\tpropNtrlM2\n"
if(model=="SC"):
	if(Nvariation=="homo"):
		if(Mvariation=="homo"):
			res="N1\tN2\tNa\tTsplit\tTsc\tM1\tM2\n"
		if(Mvariation=="hetero"):
			res="N1\tN2\tNa\tTsplit\tTsc\tM1\tM2\tshape1M1\tshape2M1\tshape1M2\tshape2M2\tpropNtrlM1\tpropNtrlM2\n"
	if(Nvariation=="hetero"):
		if(Mvariation=="homo"):
			res="N1\tN2\tNa\tTsplit\tTsc\tshape1Ne\tshape2Ne\tpropNtrlNe1\tpropNtrlNe2\tM1\tM2\n"
		if(Mvariation=="hetero"):
			res="N1\tN2\tNa\tTsplit\tTsc\tshape1Ne\tshape2Ne\tpropNtrlNe1\tpropNtrlNe2\tM1\tM2\tshape1M1\tshape2M1\tshape1M2\tshape2M2\tpropNtrlM1\tpropNtrlM2\n"

for i in range(nreps):
	n1prior=uniform(n1[0], n1[1])
	n2prior=uniform(n2[0], n2[1])
	nAprior=uniform(nA[0], nA[1])
	Tsplit=uniform(tau[0], tau[1])
	Tsmall=uniform(min(tau), Tsplit)
	M1prior=uniform(M1[0], M1[1])
	M2prior=uniform(M2[0], M2[1])
	shape1mig1=uniform(shape1[0], shape1[1])
	shape2mig1=uniform(shape2[0], shape2[1])
	shape1mig2=uniform(shape1[0], shape1[1])
	shape2mig2=uniform(shape2[0], shape2[1])
	shape1Ne=uniform(shape1[0], shape1[1])
	shape2Ne=uniform(shape2[0], shape2[1])
	TsplitGenomic=[Tsplit]*nlocus
	TsmallGenomic=[Tsmall]*nlocus
	if(Nvariation=="hetero"):
		n1priorGenomic=binomBeta(nlocus=nlocus, shape1=shape1Ne, shape2=shape2Ne, scalar=n1prior)
		n2priorGenomic=binomBeta(nlocus=nlocus, shape1=shape1Ne, shape2=shape2Ne, scalar=n2prior)
		nApriorGenomic=binomBeta(nlocus=nlocus, shape1=shape1Ne, shape2=shape2Ne, scalar=nAprior)
	if(Nvariation=="homo"):
		n1priorGenomic={}
		n2priorGenomic={}
		nApriorGenomic={}
		n1priorGenomic["values"]=[n1prior]*nlocus
		n2priorGenomic["values"]=[n2prior]*nlocus
		nApriorGenomic["values"]=[nAprior]*nlocus
	if(Mvariation=="hetero"):
		M1priorGenomic=binomBeta(nlocus=nlocus, shape1=shape1mig1, shape2=shape2mig1, scalar=M1prior)
		if(sym=="sym"):
			M2priorGenomic={}
			M2priorGenomic=M1priorGenomic["values"]
			M2prior=M1prior
		if(sym=="asym"):
			M2priorGenomic=binomBeta(nlocus=nlocus, shape1=shape1mig2, shape2=shape2mig2, scalar=M2prior)
	if(Mvariation=="homo"):
		M1priorGenomic={}
		M1priorGenomic["values"]=[M1prior]*nlocus
		if(sym=="sym"):
			M2priorGenomic={}
			M2priorGenomic["values"]=M1priorGenomic["values"]
			M2prior=M1prior
		if(sym=="asym"):
			M2priorGenomic={}
			M2priorGenomic["values"]=[M2prior]*nlocus
	res+="{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}".format(n1prior, n2prior, nAprior, Tsplit)
	if(model=="SI"):
		if(Nvariation=="homo"):
			res+=" \n"
		if(Nvariation=="hetero"):
			res+=" \t{0:.5f}\t{1:.5f}\t{2}\t{3}\n".format(shape1Ne, shape2Ne, n1priorGenomic["nNeutre"], n2priorGenomic["nNeutre"])
	if(model=="IM"):
		if(Nvariation=="homo"):
			if(Mvariation=="homo"):
				res+=" \t{0:.5f}\t{1:.5f}\n".format(M1prior, M2prior)
			if(Mvariation=="hetero"):
				res+=" \t{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6}\t{7}\n".format(M1prior, M2prior, shape1mig1, shape2mig1, shape1mig2, shape2mig2, M1priorGenomic["nNeutre"], M2priorGenomic["nNeutre"])
		if(Nvariation=="hetero"):
			if(Mvariation=="homo"):
				res+=" \t{0:.5f}\t{1:.5f}\t{2}\t{3}\t{4:.5f}\t{5:.5f}\n".format(shape1Ne, shape2Ne, n1priorGenomic["nNeutre"], n2priorGenomic["nNeutre"], M1prior, M2prior)
			if(Mvariation=="hetero"):
				res+=" \t{0:.5f}\t{1:.5f}\t{2}\t{3}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10}\t{11}\n".format(shape1Ne, shape2Ne, n1priorGenomic["nNeutre"], n2priorGenomic["nNeutre"], M1prior, M2prior, shape1mig1, shape2mig1, shape1mig2, shape2mig2, M1priorGenomic["nNeutre"], M2priorGenomic["nNeutre"])
	if(model=="AM"):
		if(Nvariation=="homo"):
			if(Mvariation=="homo"):
				res+=" \t{2:.5f}\t{0:.5f}\t{1:.5f}\n".format(M1prior, M2prior, Tsmall)
			if(Mvariation=="hetero"):
				res+=" \t{8:.5f}\t{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6}\t{7}\n".format(M1prior, M2prior, shape1mig1, shape2mig1, shape1mig2, shape2mig2, M1priorGenomic["nNeutre"], M2priorGenomic["nNeutre"], Tsmall)
		if(Nvariation=="hetero"):
			if(Mvariation=="homo"):
				res+=" \t{6:.5f}\t{0:.5f}\t{1:.5f}\t{2}\t{3}\t{4:.5f}\t{5:.5f}\n".format(shape1Ne, shape2Ne, n1priorGenomic["nNeutre"], n2priorGenomic["nNeutre"], M1prior, M2prior, Tsmall)
			if(Mvariation=="hetero"):
				res+=" \t{12:.5f}\t{0:.5f}\t{1:.5f}\t{2}\t{3}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10}\t{11}\n".format(shape1Ne, shape2Ne, n1priorGenomic["nNeutre"], n2priorGenomic["nNeutre"], M1prior, M2prior, shape1mig1, shape2mig1, shape1mig2, shape2mig2, M1priorGenomic["nNeutre"], M2priorGenomic["nNeutre"], Tsmall)
	if(model=="SC"):
		if(Nvariation=="homo"):
			if(Mvariation=="homo"):
				res+=" \t{2:.5f}\t{0:.5f}\t{1:.5f}\n".format(M1prior, M2prior, Tsmall)
			if(Mvariation=="hetero"):
				res+=" \t{8:.5f}\t{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6}\t{7}\n".format(M1prior, M2prior, shape1mig1, shape2mig1, shape1mig2, shape2mig2, M1priorGenomic["nNeutre"], M2priorGenomic["nNeutre"], Tsmall)
		if(Nvariation=="hetero"):
			if(Mvariation=="homo"):
				res+=" \t{6:.5f}\t{0:.5f}\t{1:.5f}\t{2}\t{3}\t{4:.5f}\t{5:.5f}\n".format(shape1Ne, shape2Ne, n1priorGenomic["nNeutre"], n2priorGenomic["nNeutre"], M1prior, M2prior, Tsmall)
			if(Mvariation=="hetero"):
				res+=" \t{12:.5f}\t{0:.5f}\t{1:.5f}\t{2}\t{3}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10}\t{11}\n".format(shape1Ne, shape2Ne, n1priorGenomic["nNeutre"], n2priorGenomic["nNeutre"], M1prior, M2prior, shape1mig1, shape2mig1, shape1mig2, shape2mig2, M1priorGenomic["nNeutre"], M2priorGenomic["nNeutre"], Tsmall)
	for loc in range(nlocus):
		cout=""
		if(model=="SI"):
			#msnsam tbs 200 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -ej tbs 2 1 -eN tbs tbs # for 'SI'
			#cout+="{0} {1:.5f} {2:.5f} {3} {4} {5} {6} {7} {8:.5f} {9:.5f} {10:.5f} {11:.5f} {12:.5f}".format(int(nspA[loc])+int(nspB[loc]), float(theta[loc]), float(rho[loc]), int(L[loc]), int(nspA[loc]), int(nspB[loc]), 0, 0, float(n1priorGenomic["values"][loc]), float(n2priorGenomic["values"][loc]), float(TsplitGenomic[loc]), float(TsplitGenomic[loc]), float(nApriorGenomic["values"][loc]))
			#
			#msnsam tbs 200 -s 1 -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -ej tbs 2 1 -eN tbs tbs # for 'SI'
			cout+="{0} {1} {2} {3} {4} {5:.5f} {6:.5f} {7:.5f} {8:.5f} {9:.5f}".format(int(nspA[loc])+int(nspB[loc]), int(nspA[loc]), int(nspB[loc]), 0, 0, float(n1priorGenomic["values"][loc]), float(n2priorGenomic["values"][loc]), float(TsplitGenomic[loc]), float(TsplitGenomic[loc]), float(nApriorGenomic["values"][loc]))
		if(model=="IM"):
			#msnsam tbs 200 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -ej tbs 2 1 -eN tbs tbs # for 'IM'
                        #cout+="{0} {1:.5f} {2:.5f} {3} {4} {5} {6:.5f} {7:.5f} {8:.5f} {9:.5f} {10:.5f} {11:.5f} {12:.5f}".format(int(nspA[loc])+int(nspB[loc]), float(theta[loc]), float(rho[loc]), int(L[loc]), int(nspA[loc]), int(nspB[loc]), float(M1priorGenomic["values"][loc]), float(M2priorGenomic["values"][loc]), float(n1priorGenomic["values"][loc]), float(n2priorGenomic["values"][loc]), float(TsplitGenomic[loc]), float(TsplitGenomic[loc]), float(nApriorGenomic["values"][loc]))
			#
			#msnsam tbs 200 -s 1 -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -ej tbs 2 1 -eN tbs tbs # for 'IM'
                        cout+="{0} {1} {2} {3:.5f} {4:.5f} {5:.5f} {6:.5f} {7:.5f} {8:.5f} {9:.5f}".format(int(nspA[loc])+int(nspB[loc]), int(nspA[loc]), int(nspB[loc]), float(M1priorGenomic["values"][loc]), float(M2priorGenomic["values"][loc]), float(n1priorGenomic["values"][loc]), float(n2priorGenomic["values"][loc]), float(TsplitGenomic[loc]), float(TsplitGenomic[loc]), float(nApriorGenomic["values"][loc]))
		if(model=="AM"):
			#msnsam tbs 200 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 0 -m 2 1 0 -n 1 tbs -n 2 tbs -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs # for "AM"
                        #cout+="{0} {1:.5f} {2:.5f} {3} {4} {5} {6:.5f} {7:.5f} {8:.5f} {9:.5f} {13:.5f} {14:.5f} {15:.5f} {10:.5f} {11:.5f} {12:.5f}".format(int(nspA[loc])+int(nspB[loc]), float(theta[loc]), float(rho[loc]), int(L[loc]), int(nspA[loc]), int(nspB[loc]), 0, 0, float(n1priorGenomic["values"][loc]), float(n2priorGenomic["values"][loc]), float(TsplitGenomic[loc]), float(TsplitGenomic[loc]), float(nApriorGenomic["values"][loc]), TsmallGenomic[loc], float(M1priorGenomic["values"][loc]), float(M2priorGenomic["values"][loc]))
			#msnsam tbs 200 -s 1 -I 2 tbs tbs 0 -m 1 2 0 -m 2 1 0 -n 1 tbs -n 2 tbs -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs # for "AM"
                        cout+="{0} {1} {2} {3} {4} {5:.5f} {6:.5f} {10:.5f} {11:.5f} {12:.5f} {7:.5f} {8:.5f} {9:.5f}".format(int(nspA[loc])+int(nspB[loc]), int(nspA[loc]), int(nspB[loc]), 0, 0, float(n1priorGenomic["values"][loc]), float(n2priorGenomic["values"][loc]), float(TsplitGenomic[loc]), float(TsplitGenomic[loc]), float(nApriorGenomic["values"][loc]), TsmallGenomic[loc], float(M1priorGenomic["values"][loc]), float(M2priorGenomic["values"][loc]))
		if(model=="SC"):
			#msnsam tbs 200 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs # for 'SC'
                        #cout+="{0} {1:.5f} {2:.5f} {3} {4} {5} {6:.5f} {7:.5f} {8:.5f} {9:.5f} {13:.5f} {10:.5f} {11:.5f} {12:.5f}".format(int(nspA[loc])+int(nspB[loc]), float(theta[loc]), float(rho[loc]), int(L[loc]), int(nspA[loc]), int(nspB[loc]), float(M1priorGenomic["values"][loc]), float(M2priorGenomic["values"][loc]), float(n1priorGenomic["values"][loc]), float(n2priorGenomic["values"][loc]), float(TsplitGenomic[loc]), float(TsplitGenomic[loc]), float(nApriorGenomic["values"][loc]), float(TsmallGenomic[loc]))
			#msnsam tbs 200 -s 1 -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs # for 'SC'
                        cout+="{0} {1} {2} {3:.5f} {4:.5f} {5:.5f} {6:.5f} {10:.5f} {7:.5f} {8:.5f} {9:.5f}".format(int(nspA[loc])+int(nspB[loc]), int(nspA[loc]), int(nspB[loc]), float(M1priorGenomic["values"][loc]), float(M2priorGenomic["values"][loc]), float(n1priorGenomic["values"][loc]), float(n2priorGenomic["values"][loc]), float(TsplitGenomic[loc]), float(TsplitGenomic[loc]), float(nApriorGenomic["values"][loc]), float(TsmallGenomic[loc]))
		print(cout)
outputfile=open(outputParameters, "w")
outputfile.write(res)
outputfile.close()

