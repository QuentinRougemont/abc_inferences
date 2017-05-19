
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "spmain.h"
#include "nr.h"
#include "sprun.h"
/*#include <direct.h>*/

/*global variables*/
long glob_tos_node_stack;/*global variable with index of next open stack location*/
long glob_max_node_stack;/*global variable with maximum number of  nodes in the stack*/
struct type_all **glob_node_stack;	/*array of pointers to nodes in the stack*/
long glob_min_tos_node_stack;/*glob var with minimum value taken by tos*/

/*****************************************************************************/

void main(int argc,char *argv[])
{
	/*parameters of the model, taken from inputfile*/
	int *nseqAl,*nseqBl;			/*vectors of number of sequences in species A and B for locus 0 to nloc-1*/
	int nloc;			/*number of loci*/
	int *nsl;			/*vector of number of sites for locus 0 to nloc-1*/
	long ndatasets;			/*number of replicate datest in inputfile*/
	char datafilename[128];		/*generic name for the dataset files*/

	char inputfilename[128];
	char outputfilename[128];
	char datasetfile[128];			/*name of the current datasetfile*/
	char **ar;		/*vector of arguments to Msnsam*/
	int nar;		/*number of arguments to Msnsam*/
	char *argtemp;
	FILE *finputp,*fout1p,*fout2p,*fout3p,*finp;
	/*char fout1name[_MAX_PATH],fout2name[_MAX_PATH],fout3name[_MAX_PATH];
	char home_dir[_MAX_PATH];
	char replicate_dir[_MAX_PATH];*/
	int nseqAmax;				/*maximum number of sequences in species A = max length of DNA sequence over loci*/
	int nseqBmax;				/*maximum number of sequences in species A = max length of DNA sequence over loci*/
	int nsmax;				/*maximum number of sites = max length of DNA sequence over loci*/
	char ***seqAlhs;		/*Array with DNA sequences for species A, locus l, haplotype h, sites 0 to nls[l]-1*/
	char ***seqBlhs;		/*Array with DNA sequences for species B, locus l, haplotype h, sites 0 to nls[l]-1*/
	char **seqOls;		/*Array with DNA sequences for Outgroup, locus l, sites 0 to nls[l]-1*/
	struct result_poly *result_polyl;	/*vector of structure with results of classification of sites*/
	int Fst_missing;			/*number of loci for which Fst cannot be computed because no polymorphic sites*/
	int *nspolyl;		/*vector of number of polymorphic sites for locus l from 0 to nloc-1*/


	char str[SMAX], *end;	/*arguments to fgets and strtol and strtod*/
	char sspecies[SMAX];	/*name of species*/
wchar_t wstr[SMAX];	/*arguments to fgetws and strtol and strtod*/
	int iget;
	long countdatasets;		/*loop counter*/
	int countseq;			/*loop counter*/
	int countcodons;
	char sdummy[SMAX];	/*argument to sscanf*/
	int nseqtot_temp,seqlength_temp;
	int ispecies;
	char smess2[SMAX];



int seqlength;			/*length of sequences*/
int ncodons;			/*length of sequences in codons*/
int nseqtot;
int nseqspecies1,nseqspecies2,nseqspecies3;

	/*general variables*/
	int countsimul;         	/*counter of the number of simulations already done*/
	long imain,jmain,kmain,lmain;			/*counters*/
	long seed,seedsimul;	/*random number generator seed*/
	char smess[SMAX];   /*strings used to process messages */
	float mean,stdev;
	int model_type;		/*set to DOMINANT if there are dominance classes and to NORMAL otherwise*/


	/*miscellaneous*/
	FILE *fp;
	float ft1,ft2,ft3,ft4,ft5,ft6,ft7,ft8,ft9,ft10,ft11,ft12,ft13,ft14,ft20,ft21,ft22,ft23; /*for writing to outputfile */

/***********************************************************************************/

/*memory allocations and other initializations*/
	nseqAl=alloc_int_vector(MAX_NLOC);   /*vector with number of haplotypes from species A per locus from 0 to nloc-1*/
	nseqBl=alloc_int_vector(MAX_NLOC);   /*vector with number of haplotypes from species A per locus from 0 to nloc-1*/
	nsl=alloc_int_vector(MAX_NLOC);   /*vector with number of sites per locus from 0 to nloc-1*/
	ar=charmatrix(0,NARMAX,0,SMAX);	/*vector of arguments to msnsam*/

	/*	allele1=alloc_allele(N2+1);
	allele2=alloc_allele(N2+1);
	alleledeme1=alloc_matrix_allele(Ndemes,N2+1);
	alleledeme2=alloc_matrix_allele(Ndemes,N2+1);
	na_deme1=alloc_longint_vector(Ndemes);*/

/*starts initializations by getting initial conditions*/
	if(argc == 3) {					 /*the program can take inputfilename and outputfilename as arguments*/
		strcpy(inputfilename,argv[1]); 
		strcpy(outputfilename,argv[2]);
	} else { 
		strcpy(inputfilename,"spinput.txt");	/*otherwise will use these two files*/ 
		strcpy(outputfilename,"spoutput.txt");
	}
	get_initial_conditions(inputfilename,nseqAl,nseqBl,&nloc,nsl,&ndatasets,datafilename);
	print_initial_conditions(inputfilename,nseqAl,nseqBl,nloc,nsl,ndatasets,datafilename); 
	write_initial_conditions(outputfilename,inputfilename,nseqAl,nseqBl,nloc,nsl,ndatasets,datafilename); 
	write_initial_conditions(ERRORFILE,inputfilename,nseqAl,nseqBl,nloc,nsl,ndatasets,datafilename);
	/*initialize_write_polyl(outputfilename,nloc);*/
	initialize_write_ABCstat("ABCstat.txt");
	/*initialize_write_ABCstat("ABC_param.txt");*/


	/*memory allocations and other initializations*/
	nseqAmax=nseqAl[0];
	nseqBmax=nseqBl[0];
	nsmax=nsl[0];
	for(imain=1;imain<nloc;imain++) {	/* loop over loci to search for max number sites within a sequence*/
		if(nseqAl[imain]>nseqAmax) nseqAmax=nseqAl[imain];
		if(nseqBl[imain]>nseqBmax) nseqBmax=nseqBl[imain];
		if(nsl[imain]>nsmax) nsmax=nsl[imain];
	}
	seqAlhs=char3tensor(0,nloc-1,0,nseqAmax,0,nsmax+1);
	seqBlhs=char3tensor(0,nloc-1,0,nseqBmax,0,nsmax+1);
	seqOls=charmatrix(0,nloc-1,0,nsmax+1);
	result_polyl=alloc_result_poly(nloc+2);
	nspolyl=alloc_int_vector(nloc);   /*vector with number of polymorphic sites per locus from 0 to nloc-1*/


	/* Get the current working directory: */
	/*if( _getcwd( home_dir, _MAX_PATH ) == NULL ) {
		sprintf(smess,"error in getting the home directory");
		write(ERRORFILE,smess);
		exit(1);
	}*/
	if((finputp=fopen(datafilename,"r"))!=NULL)
	{
		/*starts by reading the command line of the MS run and write it to ABC_param.txt*/
	/*	if(fgets(str,SMAX-2,finputp)) write_onlytofile("ABC_param.txt",str);
		else {
			sprintf(smess,"error in reading dataset : cannot read datasetfile ");
			write(ERRORFILE,smess);
			exit(1);
		}*/

	/*Starts loop over each replicate dataset*/		
		for(countdatasets=0;countdatasets<ndatasets;++countdatasets) {    /*loop replicate datasets*/
			get_time_short(smess2);	 
			sprintf(smess,"\nReplicate: %4ld Time: %s",countdatasets,smess2);
			imain=ndatasets/100;
			if(imain>0)	if((countdatasets % imain)==0) write(ERRORFILE,smess);
			finp=finputp;	
			result_polyl[nloc].totsites=0;	/*for sums of stat*/
			result_polyl[nloc].bialsites=0;
			result_polyl[nloc].multisites=0;
			result_polyl[nloc].sf=0;
//			result_polyl[nloc].sfB=0;
			result_polyl[nloc].sfout=0;
			result_polyl[nloc].sxA=0;
			result_polyl[nloc].sxB=0;
//			result_polyl[nloc].sxAfB=0;
//			result_polyl[nloc].sxBfA=0;
			result_polyl[nloc].ss=0;
			result_polyl[nloc].Wald=0;
			result_polyl[nloc].piA=0.0F;
			result_polyl[nloc].piB=0.0F;
			result_polyl[nloc].thetaA=0.0F;
			result_polyl[nloc].thetaB=0.0F;
			result_polyl[nloc].DA=0.0F;
			result_polyl[nloc].DB=0.0F;
			result_polyl[nloc].dAB=0.0F;
			result_polyl[nloc].dnAB=0.0F;
			result_polyl[nloc].FST=0.0F;
			result_polyl[nloc+1].totsites=0;	/*for sums of squares of stat*/
			result_polyl[nloc+1].bialsites=0;
			result_polyl[nloc+1].multisites=0;
			result_polyl[nloc+1].sf=0;
//			result_polyl[nloc+1].sfB=0;
			result_polyl[nloc+1].sfout=0;
			result_polyl[nloc+1].sxA=0;
			result_polyl[nloc+1].sxB=0;
//			result_polyl[nloc+1].sxAfB=0;
//			result_polyl[nloc+1].sxBfA=0;
			result_polyl[nloc+1].ss=0;
			result_polyl[nloc+1].Wald=0;
			result_polyl[nloc+1].piA=0.0F;
			result_polyl[nloc+1].piB=0.0F;
			result_polyl[nloc+1].thetaA=0.0F;
			result_polyl[nloc+1].thetaB=0.0F;
			result_polyl[nloc+1].DA=0.0F;
			result_polyl[nloc+1].DB=0.0F;
			result_polyl[nloc+1].dAB=0.0F;
			result_polyl[nloc+1].dnAB=0.0F;
			result_polyl[nloc+1].FST=0.0F;
			/*sprintf(datasetfile,"%s_%d.arp",datafilename,countdatasets);*/
			get_dataset(finp,countdatasets,nseqAl,nseqBl,nloc,nsl,nsmax,seqAlhs,seqBlhs,seqOls,nspolyl);
			Fst_missing=0;
			for(lmain=0;lmain<nloc;lmain++)	{ /*loop over loci */
	/*write_dataset_fasta("data_fasta.fsa",nseqAl[lmain],nseqBl[lmain],lmain,nsl[lmain],seqAlhs,seqBlhs,seqOls);*/
				compute_polyl(nseqAl[lmain],nseqBl[lmain],nsl[lmain],nspolyl[lmain],seqAlhs[lmain],seqBlhs[lmain],seqOls[lmain],&result_polyl[lmain]);
				result_polyl[nloc].totsites+=result_polyl[lmain].totsites;
				result_polyl[nloc].bialsites+=result_polyl[lmain].bialsites;
				result_polyl[nloc].multisites+=result_polyl[lmain].multisites;
				result_polyl[nloc].sf+=result_polyl[lmain].sf;
//				result_polyl[nloc].sfB+=result_polyl[lmain].sfB;
				result_polyl[nloc].sfout+=result_polyl[lmain].sfout;
				result_polyl[nloc].sxA+=result_polyl[lmain].sxA;
				result_polyl[nloc].sxB+=result_polyl[lmain].sxB;
//				result_polyl[nloc].sxAfB+=result_polyl[lmain].sxAfB;
//				result_polyl[nloc].sxBfA+=result_polyl[lmain].sxBfA;
				result_polyl[nloc].ss+=result_polyl[lmain].ss;
				result_polyl[nloc].Wald+=result_polyl[lmain].Wald;
				result_polyl[nloc].piA+=result_polyl[lmain].piA;
				result_polyl[nloc].piB+=result_polyl[lmain].piB;
				result_polyl[nloc].thetaA+=result_polyl[lmain].thetaA;
				result_polyl[nloc].thetaB+=result_polyl[lmain].thetaB;
				result_polyl[nloc].DA+=result_polyl[lmain].DA;
				result_polyl[nloc].DB+=result_polyl[lmain].DB;
				result_polyl[nloc].dAB+=result_polyl[lmain].dAB;
				result_polyl[nloc].dnAB+=result_polyl[lmain].dnAB;
				if(result_polyl[lmain].FST != MISSING) result_polyl[nloc].FST+=result_polyl[lmain].FST;
					else Fst_missing++;
				result_polyl[nloc+1].totsites+=sqr(result_polyl[lmain].totsites);
				result_polyl[nloc+1].bialsites+=sqr(result_polyl[lmain].bialsites);
				result_polyl[nloc+1].multisites+=sqr(result_polyl[lmain].multisites);
				result_polyl[nloc+1].sf+=sqr(result_polyl[lmain].sf);
//				result_polyl[nloc+1].sfB+=sqr(result_polyl[lmain].sfB);
				result_polyl[nloc+1].sfout+=sqr(result_polyl[lmain].sfout);
				result_polyl[nloc+1].sxA+=sqr(result_polyl[lmain].sxA);
				result_polyl[nloc+1].sxB+=sqr(result_polyl[lmain].sxB);
//				result_polyl[nloc+1].sxAfB+=sqr(result_polyl[lmain].sxAfB);
//				result_polyl[nloc+1].sxBfA+=sqr(result_polyl[lmain].sxBfA);
				result_polyl[nloc+1].ss+=sqr(result_polyl[lmain].ss);
				result_polyl[nloc+1].Wald+=sqr(result_polyl[lmain].Wald);
				result_polyl[nloc+1].piA+=sqr(result_polyl[lmain].piA);
				result_polyl[nloc+1].piB+=sqr(result_polyl[lmain].piB);
				result_polyl[nloc+1].thetaA+=sqr(result_polyl[lmain].thetaA);
				result_polyl[nloc+1].thetaB+=sqr(result_polyl[lmain].thetaB);
				result_polyl[nloc+1].DA+=sqr(result_polyl[lmain].DA);
				result_polyl[nloc+1].DB+=sqr(result_polyl[lmain].DB);
				result_polyl[nloc+1].dAB+=sqr(result_polyl[lmain].dAB);
				result_polyl[nloc+1].dnAB+=sqr(result_polyl[lmain].dnAB);
				if(result_polyl[lmain].FST != MISSING) result_polyl[nloc+1].FST+=sqr(result_polyl[lmain].FST);
			}	/*end loop over loci*/
			result_polyl[nloc+1].totsites=comp_std(nloc,result_polyl[nloc].totsites,result_polyl[nloc+1].totsites);
			result_polyl[nloc+1].bialsites=comp_std(nloc,result_polyl[nloc].bialsites,result_polyl[nloc+1].bialsites);
			result_polyl[nloc+1].multisites=comp_std(nloc,result_polyl[nloc].multisites,result_polyl[nloc+1].multisites);
			result_polyl[nloc+1].sf=comp_std(nloc,result_polyl[nloc].sf,result_polyl[nloc+1].sf);
//			result_polyl[nloc+1].sfB=comp_std(nloc,result_polyl[nloc].sfB,result_polyl[nloc+1].sfB);
			result_polyl[nloc+1].sfout=comp_std(nloc,result_polyl[nloc].sfout,result_polyl[nloc+1].sfout);
			result_polyl[nloc+1].sxA=comp_std(nloc,result_polyl[nloc].sxA,result_polyl[nloc+1].sxA);
			result_polyl[nloc+1].sxB=comp_std(nloc,result_polyl[nloc].sxB,result_polyl[nloc+1].sxB);
//			result_polyl[nloc+1].sxAfB=comp_std(nloc,result_polyl[nloc].sxAfB,result_polyl[nloc+1].sxAfB);
//			result_polyl[nloc+1].sxBfA=comp_std(nloc,result_polyl[nloc].sxBfA,result_polyl[nloc+1].sxBfA);
			result_polyl[nloc+1].ss=comp_std(nloc,result_polyl[nloc].ss,result_polyl[nloc+1].ss);
			result_polyl[nloc+1].Wald=comp_std(nloc,result_polyl[nloc].Wald,result_polyl[nloc+1].Wald);
			result_polyl[nloc+1].piA=comp_std(nloc,result_polyl[nloc].piA,result_polyl[nloc+1].piA);
			result_polyl[nloc+1].piB=comp_std(nloc,result_polyl[nloc].piB,result_polyl[nloc+1].piB);
			result_polyl[nloc+1].thetaA=comp_std(nloc,result_polyl[nloc].thetaA,result_polyl[nloc+1].thetaA);
			result_polyl[nloc+1].thetaB=comp_std(nloc,result_polyl[nloc].thetaB,result_polyl[nloc+1].thetaB);
			result_polyl[nloc+1].DA=comp_std(nloc,result_polyl[nloc].DA,result_polyl[nloc+1].DA);
			result_polyl[nloc+1].DB=comp_std(nloc,result_polyl[nloc].DB,result_polyl[nloc+1].DB);
			result_polyl[nloc+1].dAB=comp_std(nloc,result_polyl[nloc].dAB,result_polyl[nloc+1].dAB);
			result_polyl[nloc+1].dnAB=comp_std(nloc,result_polyl[nloc].dnAB,result_polyl[nloc+1].dnAB);
			if((nloc-Fst_missing)>0) result_polyl[nloc+1].FST=comp_std((nloc-Fst_missing),result_polyl[nloc].FST,result_polyl[nloc+1].FST);
			else result_polyl[nloc+1].FST=0;
			result_polyl[nloc].totsites /= (float) nloc;
			result_polyl[nloc].bialsites /= (float) nloc;
			result_polyl[nloc].multisites /= (float) nloc;
			result_polyl[nloc].sf /= (float) nloc;
//			result_polyl[nloc].sfB /= (float) nloc;
			result_polyl[nloc].sfout /= (float) nloc;
			result_polyl[nloc].sxA /= (float) nloc;
			result_polyl[nloc].sxB /= (float) nloc;
//			result_polyl[nloc].sxAfB /= (float) nloc;
//			result_polyl[nloc].sxBfA /= (float) nloc;
			result_polyl[nloc].ss /= (float) nloc;
			result_polyl[nloc].Wald /= (float) nloc;
			result_polyl[nloc].piA /= (float) nloc;
			result_polyl[nloc].piB /= (float) nloc;
			result_polyl[nloc].thetaA /= (float) nloc;
			result_polyl[nloc].thetaB /= (float) nloc;
			result_polyl[nloc].DA /= (float) nloc;
			result_polyl[nloc].DB /= (float) nloc;
			result_polyl[nloc].dAB /= (float) nloc;
			result_polyl[nloc].dnAB /= (float) nloc;
			if((nloc-Fst_missing)>0) result_polyl[nloc].FST /= (float) (nloc-Fst_missing);
			else result_polyl[nloc].FST=0;
			/*write_polyl(outputfilename,countdatasets,nloc,result_polyl);*/
			write_ABCstat("ABCstat.txt",countdatasets,nloc,result_polyl);
		}	/*end loop over replicate datasets*/
	fclose(finputp);
	}     /*if can open filename*/
		else {    /*if cannot open filename*/
			sprintf(smess,"\nError in reading data set file :\n\tCannot open file %s",datafilename);
			write(ERRORFILE,smess);
			exit(1);
	}/*end of else*/
	exit(1);
}    /*end Main*/
