#include <stdio.h>
#include <stdlib.h>
#include "spmain.h"
#include <math.h>
#include <string.h>

void get_initial_conditions(const char *inputfile,int *nseqAl,int *nseqBl,int *nloc, int *nsl, long *ndatasets, char *datafile)
 {
/* reads inputfile for values of parameters and check them against legal range*/
	FILE *fp;
	char str[SMAX], *end;	/*arguments to fgets and strtol and strtod*/
	char swhiteline[SMAX];
	int l;
	int iget;
	long lget;
	double dget;
	char smess[SMAX];

	if((fp=fopen(inputfile,"r"))!=NULL)
	{
		do {
			if(fgets(str,SMAX-2,fp));
			else {
				sprintf(smess,"error in reading initial conditions : cannot read datafile");
				write(ERRORFILE,smess);
				exit(1);
			}
			iget=sscanf(str,"%s",swhiteline);
		} while(iget<1);
		iget=sscanf(str,"%d",nloc);
		if(iget<1) {                              /*test nloc*/
			sprintf(smess,"error in reading initial conditions : cannot read nloc"); 
			write(ERRORFILE,smess);
			exit(1);
		}
		for(l=0;l<*nloc;l++) {	/*loop over loci*/
			if(fgets(str,SMAX-2,fp));	/*read nseqAl at locus l*/
			else {
				sprintf(smess,"error in reading initial conditions : cannot read nseqAl at locus %d",l);
				write(ERRORFILE,smess);
				exit(1);
			}
			iget=sscanf(str,"%d",&nseqAl[l]);		/*read nseqA*/
			if(iget<1) {                              /*test nseqA*/
				sprintf(smess,"error in reading initial conditions : cannot read nseqAl at locus %d",l); 
				write(ERRORFILE,smess);
				exit(1);
			} 
			if(fgets(str,SMAX-2,fp));	/*read nseqB*/
			else {
				sprintf(smess,"error in reading initial conditions : cannot read nseqBl at locus %d",l);
				write(ERRORFILE,smess);
				exit(1);
			}
			iget=sscanf(str,"%d",&nseqBl[l]);
			if(iget<1) {                              /*test nseqB*/
				sprintf(smess,"error in reading initial conditions : cannot read nseqBl at locus %d",l); 
				write(ERRORFILE,smess);
				exit(1);
			} 

			if(fgets(str,SMAX-2,fp));	/*read nsl number of sites at locus l*/
			else {
				sprintf(smess,"error in reading initial conditions : cannot read numb sites at locus %d",l);
				write(ERRORFILE,smess);
				exit(1);
			}
			iget=sscanf(str,"%d",&nsl[l]);
			if(iget<1) {                              /*test nsl number of sites at locus l*/
				sprintf(smess,"error in reading initial conditions : cannot read numb sites at locus %d",l); 
				write(ERRORFILE,smess);
				exit(1);
			}
		}	/*end loop over loci*/
		if(fgets(str,SMAX-2,fp));	/*read ndatasets*/
		else {
			sprintf(smess,"error in reading initial conditions : cannot read ndatasets");
			write(ERRORFILE,smess);
			exit(1);
		}
		iget=sscanf(str,"%ld",ndatasets);
		if(iget<1) {                              /*test ndatasets*/
			sprintf(smess,"error in reading initial conditions : cannot read ndatasets"); 
			write(ERRORFILE,smess);
			exit(1);
		}
		if(fgets(str,SMAX-2,fp));	/*read datafile*/
		else {
			sprintf(smess,"error in reading initial conditions : cannot read nseqA");
			write(ERRORFILE,smess);
			exit(1);
		}
		iget=sscanf(str,"%s",datafile);
		if(iget<1) {                              /*test datafile*/
			sprintf(smess,"error in reading initial conditions : cannot read datafile"); 
			write(ERRORFILE,smess);
			exit(1);
		}
		fclose(fp);
	}     /*if can open filename*/
		else {    /*if cannot open filename*/
			sprintf(smess,"\nError in reading initial conditions :\n\tCannot open file %s",inputfile);
			write(ERRORFILE,smess);
			exit(1);
	}
}		 /*end of get_initial_conditions*/

/*****************************************************************************/

void get_dataset(FILE *fp,long dataset,int *nseqAl,int *nseqBl,int nloc,int *nsl,int nsmax,char ***seqAlhs,char ***seqBlhs,char **seqOls, int *nspolyl)
{
/* reads datasetfile and extracts DNA sequences from species A, B and the outgroup*/
	char str[SMAX], *end;	/*arguments to fgets and strtol and strtod*/
	char swhiteline[SMAX];
	int l,h,copies,pos,nseg,nnonseg;
	int iget;
	long lget;
	double dget;
	char smess[SMAX];
	char seql[SMAX];
	char haplotype[SMAX];
	char seqanc[SMAX];
	/*int param1,param4,param5,param6;
	float param2,param3,param7, param8, param9,param10, param11,param12;
	float previous_param12;*/ 

	for(pos=0;pos<SMAX-1;pos++) seqanc[pos]='0';
	for(l=0;l<nloc;l++) {	/*loop over loci*/
		if(l==0) {	/*for the first locus, get and print to file the parameter values*/
			/*first get the parameter values*/
			do {						/*loop over lines in file until found // beginning the parameter line*/
				if(fgets(str,SMAX-2,fp));
				else {
					sprintf(smess,"error in reading dataset %ld : cannot read datasetfile ",dataset);
					write(ERRORFILE,smess);
					exit(1);
				}
				swhiteline[0]='\0';
				iget=sscanf(str,"%s",swhiteline); /*read characters until white space*/
			} while(strcmp(swhiteline,"//") != 0); /*found parameter line */
/*			iget=sscanf(str+3,"%d %f %f %d %d %d %f %f %f %f %f %f",&param1,&param2,&param3,&param4,&param5,&param6,
						&param7,&param8,&param9,&param10,&param11,&param12);	get the parameters values
			if(iget<12) {
				sprintf(smess,"error in reading dataset %ld : parameters not read successfully ",dataset);
				write(ERRORFILE,smess);
				sprintf(smess,"Dataset %ld Line: %s",dataset,str);
				write(ERRORFILE,smess);
				exit(1);
			}
			previous_param12 = param12;
			sprintf(smess,"\n%ld\t%f\t%f\t%f\t%f\t%f\t%f",dataset,param7,param8,param9,param10,param11,param12);
			write_onlytofile("ABC_param.txt",smess);*/
			/*now get the number of seg sites for locus 0*/
			if(fgets(str,SMAX-2,fp));
			else {
				sprintf(smess,"error in reading dataset %ld : cannot read datasetfile ",dataset);
				write(ERRORFILE,smess);
				exit(1);
			}
			iget=sscanf(str+10,"%d",&nseg);	/*get the number of seg sites*/
			if(iget<1) {
				sprintf(smess,"error in reading dataset %ld : number of seg sites not read successfully ",dataset);
				write(ERRORFILE,smess);
				sprintf(smess,"Dataset %ld Line: %s",dataset,str);
				write(ERRORFILE,smess);
				exit(1);
			}
		} else {	/*of of case with l==0 start of case with l >0*/
			do {						/*loop over lines in file until found // beginning the parameter line*/
				if(fgets(str,SMAX-2,fp));
				else {
					sprintf(smess,"error in reading dataset %ld : cannot read datasetfile ",dataset);
					write(ERRORFILE,smess);
					exit(1);
				}
				swhiteline[0]='\0';
				iget=sscanf(str,"%s",swhiteline);
			} while(strcmp(swhiteline,"//") != 0); /*found parameter line */
	/*		iget=sscanf(str+3,"%d %f %f %d %d %d %f %f %f %f %f %f",&param1,&param2,&param3,&param4,&param5,&param6,
						&param7,&param8,&param9,&param10,&param11,&param12);	get the parameters values
			if(iget<12) {
				sprintf(smess,"error in reading dataset : parameters not read successfully ");
				write(ERRORFILE,smess);
				sprintf(smess,"Dataset %ld Line: %s",dataset,str);
				write(ERRORFILE,smess);
				exit(1);
			}
			if(param12 != previous_param12) {
				sprintf(smess,"error in reading dataset %ld: Last Parameter locus %d = %f whereas Last Parameter locus 0 = %f",dataset,l,previous_param12,param12);
				write(ERRORFILE,smess);
				sprintf(smess,"Dataset %ld Line: %s",dataset,str);
				write(ERRORFILE,smess);
				exit(1);
			}*/
			/*now get the number of seg sites for locus 0*/
			if(fgets(str,SMAX-2,fp));
			else {
				sprintf(smess,"error in reading dataset %ld : cannot read datasetfile ",dataset);
				write(ERRORFILE,smess);
				exit(1);
			}
			iget=sscanf(str+10,"%d",&nseg);	/*get the number of seg sites*/
			if(iget<1) {
				sprintf(smess,"error in reading dataset %ld : number of seg sites not read successfully ",dataset);
				write(ERRORFILE,smess);
				sprintf(smess,"Line %s",str);
				write(ERRORFILE,smess);
				exit(1);
			}

		}	/*end of else for l>0*/
		nspolyl[l]=nseg;
		if(nseg==0) {	/*when no segregating sites are present in the dataset*/ 
			for(h=0;h<nseqAl[l];h++) {	/*loop over haplotypes in species A*/
				strncpy(seqAlhs[l][h],seqanc,nsl[l]);
				seqAlhs[l][h][nsl[l]]='\0';
			}	/* end loop over haplotypes in species A*/
			for(h=0;h<nseqBl[l];h++) {	/*loop over haplotypes in species B*/
				strncpy(seqBlhs[l][h],seqanc,nsl[l]);
				seqBlhs[l][h][nsl[l]]='\0';
			}	/* end loop over haplotypes in species B*/
		} else {	/*if at least one segregating site*/
			if(nseg<=nsl[l]) nnonseg=nsl[l]-nseg;
			else {	/*if too much segregating sites*/
				nspolyl[l]=nsl[l];
				nnonseg=0;
				sprintf(smess,"/nerror in reading dataset : nseg=%d > total number of sites=%d in locus %d",nseg,nsl[l],l);
				write_onlytofile(ERRORFILE,smess);
			}	/*of else*/
			if(fgets(str,SMAX-2,fp));	/*skip the line with seg sites positions*/
			/*now starts Species A*/
			for(h=0;h<nseqAl[l];h++) {	/*loop over haplotypes in species A*/
				if(fgets(str,SMAX-2,fp));
				else {
					sprintf(smess,"error in reading dataset : cannot read seq of haplotype %d in species A",h);
					write(ERRORFILE,smess);
					exit(1);
				}
				strcpy(seqAlhs[l][h],str);
				strncpy(seqAlhs[l][h]+nseg,seqanc,nnonseg);
				seqAlhs[l][h][nsl[l]]='\0';
			}	/* end loop over haplotypes in species A*/
			for(h=0;h<nseqBl[l];h++) {	/*loop over haplotypes in species B*/
				if(fgets(str,SMAX-2,fp));
				else {
					sprintf(smess,"error in reading dataset : cannot read seq of haplotype %d in species B",h);
					write(ERRORFILE,smess);
					exit(1);
				}
				strcpy(seqBlhs[l][h],str);
				strncpy(seqBlhs[l][h]+nseg,seqanc,nnonseg);
				seqBlhs[l][h][nsl[l]]='\0';
			}	/* end loop over haplotypes in species B*/
		}	/*end of else if at least one segregating site*/
		strncpy(seqOls[l],seqanc,nsl[l]);
		seqOls[l][nsl[l]]='\0';
	}	/*end loop over loci*/
}		 /*end of get_dataset*/

/*****************************************************************************/
void get_dataset_bug(FILE *fp,long dataset,int *nseqAl,int *nseqBl,int nloc,int *nsl,int nsmax,char ***seqAlhs,char ***seqBlhs,char **seqOls, int *nspolyl)
/* reads datasetfile and extracts DNA sequences from species A, B and the outgroup*/
{
	char str[SMAX], *end;	/*arguments to fgets and strtol and strtod*/
	char swhiteline[SMAX];
	int l,h,copies,pos,nseg,nnonseg;
	int iget;
	long lget;
	double dget;
	char smess[SMAX];
	char seql[SMAX];
	char haplotype[SMAX];
	char seqanc[SMAX];
	int param1,param4,param5,param6;
	float param2,param3,param7, param8, param9,param10, param11,param12;
	float previous_param12;

	for(pos=0;pos<SMAX-1;pos++) seqanc[pos]='0';
	for(l=0;l<nloc;l++) {	/*loop over loci*/
		if(l==0) {	/*for the first locus, get and print to file the parameter values*/
			/*first get the parameter values*/
			do {						/*loop over lines in file until found // beginning the parameter line*/
				if(fgets(str,SMAX-2,fp));
				else {
					sprintf(smess,"error in reading dataset %ld : cannot read datasetfile ",dataset);
					write(ERRORFILE,smess);
					exit(1);
				}
				swhiteline[0]='\0';
				iget=sscanf(str,"%s",swhiteline);
			} while(strcmp(swhiteline,"//") != 0); /*found parameter line */
/*			iget=sscanf(str+3,"%d %f %f %d %d %d %f %f %f %f %f %f",&param1,&param2,&param3,&param4,&param5,&param6,
						&param7,&param8,&param9,&param10,&param11,&param12);*/	/*get the parameters values*/
/*			if(iget<12) {
				sprintf(smess,"error in reading dataset %ld : parameters not read successfully ",dataset);
				write(ERRORFILE,smess);
				sprintf(smess,"Dataset %ld Line: %s",dataset,str);
				write(ERRORFILE,smess);
				exit(1);
			}
			previous_param12 = param12;
			sprintf(smess,"\n%ld\t%f\t%f\t%f\t%f\t%f\t%f",dataset,param7,param8,param9,param10,param11,param12);
			write_onlytofile("ABC_param.txt",smess);*/
			/*now get the number of seg sites for locus 0*/
			if(fgets(str,SMAX-2,fp));
			else {
				sprintf(smess,"error in reading dataset %ld : cannot read datasetfile ",dataset);
				write(ERRORFILE,smess);
				exit(1);
			}
			iget=sscanf(str+10,"%d",&nseg);	/*get the number of seg sites*/
			if(iget<1) {
				sprintf(smess,"error in reading dataset %ld : number of seg sites not read successfully ",dataset);
				write(ERRORFILE,smess);
				sprintf(smess,"Dataset %ld Line: %s",dataset,str);
				write(ERRORFILE,smess);
				exit(1);
			}
		} else {	/*of of case with l==0 start of case with l >0*/
			do {						/*loop over lines in file until found // beginning the parameter line*/
				if(fgets(str,SMAX-2,fp));
				else {
					sprintf(smess,"error in reading dataset %ld : cannot read datasetfile ",dataset);
					write(ERRORFILE,smess);
					exit(1);
				}
				swhiteline[0]='\0';
				iget=sscanf(str,"%s",swhiteline);
			} while(strcmp(swhiteline,"//") != 0); /*found parameter line */
			iget=sscanf(str+3,"%d %f %f %d %d %d %f %f %f %f %f %f",&param1,&param2,&param3,&param4,&param5,&param6,
						&param7,&param8,&param9,&param10,&param11,&param12);	/*get the parameters values*/
/*			if(iget<12) {
				sprintf(smess,"error in reading dataset : parameters not read successfully ");
				write(ERRORFILE,smess);
				sprintf(smess,"Dataset %ld Line: %s",dataset,str);
				write(ERRORFILE,smess);
				exit(1);
			}*/
/*			if(param12 != previous_param12) {
				sprintf(smess,"error in reading dataset %ld: Last Parameter locus %d = %f whereas Last Parameter locus 0 = %f",dataset,l,previous_param12,param12);
				write(ERRORFILE,smess);
				sprintf(smess,"Dataset %ld Line: %s",dataset,str);
				write(ERRORFILE,smess);
				exit(1);
			}*/
			/*now get the number of seg sites for locus 0*/
			if(fgets(str,SMAX-2,fp));
			else {
				sprintf(smess,"error in reading dataset %ld : cannot read datasetfile ",dataset);
				write(ERRORFILE,smess);
				exit(1);
			}
			iget=sscanf(str+10,"%d",&nseg);	/*get the number of seg sites*/
			if(iget<1) {
				sprintf(smess,"error in reading dataset %ld : number of seg sites not read successfully ",dataset);
				write(ERRORFILE,smess);
				sprintf(smess,"Line %s",str);
				write(ERRORFILE,smess);
				exit(1);
			}

		}	/*end of else for l>0*/
		nspolyl[l]=nseg;
		if(nseg==0) {
			for(h=0;h<nseqAl[l];h++) {	/*loop over haplotypes in species A*/
				strncpy(seqAlhs[l][h],seqanc,nsl[l]);
				seqAlhs[l][h][nsl[l]]='\0';
			}	/* end loop over haplotypes in species A*/
			for(h=0;h<nseqBl[l];h++) {	/*loop over haplotypes in species B*/
				strncpy(seqBlhs[l][h],seqanc,nsl[l]);
				seqBlhs[l][h][nsl[l]]='\0';
			}	/* end loop over haplotypes in species B*/
		} else {	/*if at least one segregating site*/
			if(nseg<=nsl[l]) nnonseg=nsl[l]-nseg;
			else {	/*if too much segregating sites*/
				nspolyl[l]=nsl[l];
				nnonseg=0;
				sprintf(smess,"/nerror in reading dataset : nseg=%d > total number of sites=%d in locus %d",nseg,nsl[l],l);
				write_onlytofile(ERRORFILE,smess);
			}	/*of else*/
			if(fgets(str,SMAX-2,fp));	/*skip one line*/
			/*now starts Species A*/
			for(h=0;h<nseqAl[l];h++) {	/*loop over haplotypes in species A*/
				if(fgets(str,SMAX-2,fp));
				else {
					sprintf(smess,"error in reading dataset : cannot read seq of haplotype %d in species A",h);
					write(ERRORFILE,smess);
					exit(1);
				}
				strcpy(seqAlhs[l][h],str);
				strncpy(seqAlhs[l][h]+nseg,seqanc,nnonseg);
				seqAlhs[l][h][nsl[l]]='\0';
			}	/* end loop over haplotypes in species A*/
			for(h=0;h<nseqBl[l];h++) {	/*loop over haplotypes in species B*/
				if(fgets(str,SMAX-2,fp));
				else {
					sprintf(smess,"error in reading dataset : cannot read seq of haplotype %d in species B",h);
					write(ERRORFILE,smess);
					exit(1);
				}
				strcpy(seqBlhs[l][h],str);
				strncpy(seqBlhs[l][h]+nseg,seqanc,nnonseg);
				seqBlhs[l][h][nsl[l]]='\0';
			}	/* end loop over haplotypes in species B*/
		}	/*end of else if at least one segregating site*/
		strncpy(seqOls[l],seqanc,nsl[l]);
		seqOls[l][nsl[l]]='\0';
	}	/*end loop over loci*/
}		 /*end of get_dataset_bug*/


/*********************************************************************************************/

void get_initial_conditions_dynamics(const char *inputfile,int *model_of_reproduction,
								long *N,int *Ndemes, double *geneflow,
								double *mutrate1,double *mutrate2,
								double *recombrate,double *ovd_s,
								int *nsimul,int *na_sample,
								long *nall_life,int *nclas_life,float *freq_life)
 {
/* reads inputfile for values of parameters and check them against legal range*/
	FILE *fp;
	char str[SMAX], *end;	/*arguments to fgets and strtol and strtod*/
	int iget,i;
	long lget;
	double dget;
	char smess[SMAX];
	if((fp=fopen(inputfile,"r"))!=NULL)
	{
		if(fgets(str,SMAX-2,fp)) iget=(int) strtol(str,&end,10);     /*read model_of_reproduction*/
			else {
				sprintf(smess,"error in reading initial conditions : cannot read model_of_reproduction");
				write(ERRORFILE,smess);
				exit(1);
			}
		if((iget<1) || (iget>MAX_VALUE_MODELS)) {                              /*test model_of_reproduction*/
			sprintf(smess,"error in reading initial conditions : model_of_reproduction=%d is invalid",iget); 
			write(ERRORFILE,smess);
			exit(1);
		} else *model_of_reproduction=iget;
		if(fgets(str,SMAX-2,fp)) lget= strtol(str,&end,10);     /*read N*/
			else {
				sprintf(smess,"error in reading initial conditions : cannot read N");
				write(ERRORFILE,smess);
				exit(1);
			}
		if((lget>NMAX) || (lget<1)){                            /*test N*/
			sprintf(smess,"error in reading initial conditions : N=%ld is higher than Nmax=%ld or invalid",lget,NMAX); 
			write(ERRORFILE,smess);
			exit(1); 
		} 	else *N=lget;
		if(fgets(str,SMAX-2,fp)) iget=(int) strtol(str,&end,10);     /*read Ndemes*/
			else {
				sprintf(smess,"error in reading initial conditions : cannot read number of demes");
				write(ERRORFILE,smess);
				exit(1);
			}
		if((iget<1) || (iget>*N)) {                              /*test Ndemes*/
			sprintf(smess,"error in reading initial conditions : Number of demes=%d is invalid",iget); 
			write(ERRORFILE,smess);
			exit(1);
		} else *Ndemes=iget;
		if(fgets(str,SMAX-2,fp)) dget=atof(str);     /*read geneflow*/
			else {
				sprintf(smess,"error in reading initial conditions : cannot read gene flow m");
				write(ERRORFILE,smess);
				exit(1);
			}
		if((dget>=1) || (dget<0)) {                     /*test geneflow*/
			sprintf(smess,"error in reading initial conditions : gene flow m=%f is invalid",dget);
			write(ERRORFILE,smess);
			exit(1);
		} else *geneflow=dget;
		if(fgets(str,SMAX-2,fp)) dget=atof(str);     /*read mutrate1*/
			else {
				sprintf(smess,"error in reading initial conditions : cannot read mutrate1");
				write(ERRORFILE,smess);
				exit(1);
			}
		if((dget>=1) || (dget<0)) {                     /*test mutrate1*/
			sprintf(smess,"error in reading initial conditions : mutrate1=%f is invalid",dget);
			write(ERRORFILE,smess);
			exit(1);
		} else *mutrate1=dget;
		if(fgets(str,SMAX-2,fp)) dget=atof(str);     /*read mutrate2*/
			else {
				sprintf(smess,"error in reading initial conditions : cannot read mutrate2");
				write(ERRORFILE,smess);
				exit(1);
			}
		if((dget>=1) || (dget<0)) {                     /*test mutrate2*/
			sprintf(smess,"error in reading initial conditions : mutrate2=%f is invalid",dget);
			write(ERRORFILE,smess);
			exit(1);
		} else *mutrate2=dget;
		if(fgets(str,126,fp)) dget=atof(str);     /*read recombrate*/
			else {
				sprintf(smess,"error in reading initial conditions : cannot read recombrate");
				exit(1);
			}
		if((dget>0.5) || (dget<0)) {                     /*test recombrate*/
			sprintf(smess,"error in reading initial conditions : recombrate=%f is invalid",dget);
			exit(1);
		} else *recombrate=dget;
		if(*model_of_reproduction==OVERDOMINANT_S) {
			if(fgets(str,126,fp)) dget=atof(str);     /*read selection coefficient*/
				else {
					sprintf(smess,"error in reading initial conditions : cannot read selection coefficient");
					write(ERRORFILE,smess);
					exit(1);
				}
			if((dget>1.0) || (dget<0)) {                     /*test selection coefficient*/
				sprintf(smess,"error in reading initial conditions : selection coefficient=%f is invalid",dget);
				write(ERRORFILE,smess);
				exit(1);
			} else *ovd_s=dget;
		}	/*end of if symmetric overdominant model*/
		if(fgets(str,SMAX-2,fp)) iget=(int) strtol(str,&end,10);     /*read nsimul*/
			else {
				sprintf(smess,"error in reading initial conditions : cannot read nsimul");
				write(ERRORFILE,smess);
				exit(1);
			}
		if((iget>NSIMULMAX) || (iget<1)) {                              /*test nsimul*/
			sprintf(smess,"error in reading initial conditions : nsimul=%d is invalid",iget); 
			write(ERRORFILE,smess);
			exit(1);
		} else *nsimul=iget;
		if(fgets(str,SMAX-2,fp)) iget=(int) strtol(str,&end,10);     /*read na_sample*/
			else {
				sprintf(smess,"error in reading initial conditions : cannot read na_sample");
				write(ERRORFILE,smess);
				exit(1);
			}
		if((iget>2*(*N)) || (iget<2)) {                              /*test na_sample*/
			sprintf(smess,"error in reading initial conditions : na_sample=%d is invalid",iget); 
			write(ERRORFILE,smess);
			exit(1);
		} else *na_sample=iget;
		if(fgets(str,SMAX-2,fp)) lget= strtol(str,&end,10);     /*read nall_life*/
			else {
				sprintf(smess,"error in reading initial conditions : cannot read nall_life");
				write(ERRORFILE,smess);
				exit(1);
			}
		if((lget>NALL_LIFE_MAX) || (lget<1)){                            /*test nall_life*/
			sprintf(smess,"error in reading initial conditions : nall_life=%ld is higher than max value=%ld or invalid",lget,NALL_LIFE_MAX); 
			write(ERRORFILE,smess);
			exit(1); 
		} 	else *nall_life=lget;
		if(fgets(str,SMAX-2,fp)) iget=(int) strtol(str,&end,10);     /*read nclas_life*/
			else {
				sprintf(smess,"error in reading initial conditions : cannot read nclas_life");
				write(ERRORFILE,smess);
				exit(1);
			}
		if((iget>NCLAS_LIFE_MAX) || (iget<1)) {                              /*test nclas_life*/
			sprintf(smess,"error in reading initial conditions : nclas_life=%d is invalid",iget); 
			write(ERRORFILE,smess);
			exit(1);
		} else *nclas_life=iget;
		for(i=0;i<*nclas_life;i++) {
			if(fgets(str,SMAX-2,fp)) {
				iget=sscanf(str,"%f",&freq_life[i]);
				if(iget != 1) {
					sprintf(smess,"error in reading initial conditions : cannot read vector of equ. freq for dynamics");
					write(ERRORFILE,smess);
					exit(1);
				}
			} else	{
					sprintf(smess,"error in reading initial conditions : cannot read vector of equ. freq for dynamics");
					write(ERRORFILE,smess);
					exit(1);
			}
		}	/*end loop over allele levels*/

		fclose(fp);
	}     /*if can open filename*/
		else {    /*if cannot open filename*/
			sprintf(smess,"\nError in reading initial conditions :\n\tCannot open file %s",inputfile);
			write(ERRORFILE,smess);
			sprintf(smess,"\nthis file should be like follows:\n\n");
			write(ERRORFILE,smess);
			sprintf(smess,"model_of_reproduction (=1 to 9)\npopulation size (<=%ld)",(long) NMAX);
			write(ERRORFILE,smess);
			sprintf(smess,"\nmutation rate\nnumber of simulations (<=%ld)",(long) NSIMULMAX); 
			write(ERRORFILE,smess);
			sprintf(smess,"\nnumber of alleles to sample\n"); 
			write(ERRORFILE,smess);
			exit(1);
	}
}		 /*end of get_initial_conditions_dynamics*/

/*****************************************************************************/

void print_initial_conditions(const char *inputfile,int *nseqAl,int *nseqBl,int nloc, int *nsl, long ndatasets, char *datafile) 
/* print initial values of parameters read 
from input file as a check*/
{
	int l;
	printf("\nThe following informations have been successfully loaded from %s:",inputfile);
	printf("\n\tnumber of loci : %d",nloc);
	for(l=0;l<nloc;l++) {	/*loop over loci*/
		printf("\n\tnumber of sequences in species A at locus %d: %d",l,nseqAl[l]);
		printf("\n\tnumber of sequences in species B at locus %d: %d",l,nseqBl[l]);
		printf("\n\t\tnumber of sites at locus %d: %d",l,nsl[l]);
	}
	printf("\n\tnumber of replicate datasets: %ld ",ndatasets); 
	printf("\n\tgeneric name of the datasets: %s",datafile);
	printf("\n\n");

}	  /*end of print_initial_conditions*/

/*****************************************************************************/

void write_initial_conditions(char *filename, const char *inputfile,int *nseqAl,int *nseqBl,int nloc, int *nsl, long ndatasets, char *datafile) 
/* write initial values of parameters read from input file to logfile as a check*/
{
	FILE *fp;
	char smess[SMAX];
	int l;

	if((fp=fopen(filename,"w"))!= NULL) {
		fprintf(fp,"The following informations have been successfully loaded from %s:",inputfile);
		fprintf(fp,"\n\tnumber of loci : %d",nloc);
		for(l=0;l<nloc;l++) {	/*loop over loci*/
		fprintf(fp,"\n\tnumber of sequences in species A at locus %d: %d",l,nseqAl[l]);
		fprintf(fp,"\n\tnumber of sequences in species B at locus %d: %d",l,nseqBl[l]);
			fprintf(fp,"\n\t\tnumber of sites at locus %d: %d",l,nsl[l]);
		}
		fprintf(fp,"\n\tnumber of replicate datasets: %ld ",ndatasets); 
		fprintf(fp,"\n\tgeneric name of the datasets: %s",datafile);
		fprintf(fp,"\n\n");
		fclose(fp);
	}     /*if can open filename*/
	else {    /*if cannot open filename*/
		sprintf(smess,"\nError in writing initial conditions :\n\tCannot open file %s",filename);
		write(ERRORFILE,smess);
		exit(1);
	}
}	  /*end of write_initial_conditions*/

/*****************************************************************************/

void write_dataset_fasta(char *filename,int nseqA,int nseqB,int loc, int nsl,char ***seqAlhs,char ***seqBlhs,char **seqOls) 
/* write a fasta file with data at locus loc  But does not stop at the 60th character as it should*/
{
	FILE *fp;
	char smess[SMAX];
	int h,s;

	if((fp=fopen(filename,"w"))!= NULL) {
		for(h=0;h<nseqA;h++) {	/*loop over haplotypes in species A*/
			fprintf(fp,"\n>SpeciesA_%d\n",h);
			for(s=0;s<nsl;s++) fprintf(fp,"%c",seqAlhs[loc][h][s]);
		}
		for(h=0;h<nseqB;h++) {	/*loop over haplotypes in species A*/
			fprintf(fp,"\n>SpeciesB_%d\n",h);
			for(s=0;s<nsl;s++) fprintf(fp,"%c",seqBlhs[loc][h][s]);
		}
		fprintf(fp,"\n>Outgroup\n");
		for(s=0;s<nsl;s++) fprintf(fp,"%c",seqOls[loc][s]);
		fprintf(fp,"\n");
		fclose(fp);
	}     /*if can open filename*/
	else {    /*if cannot open filename*/
		sprintf(smess,"\nError in  write_dataset_fasta :\n\tCannot open file %s",filename);
		write(ERRORFILE,smess);
		exit(1);
	}
}	  /*end of write_dataset_fasta*/
/*****************************************************************************/

void write_polyl(char *filename,long dataset,int nloc,struct result_poly *r)
/*write statistics of shared and fixed polymorphism to outputfile one line for each replicate dataset*/
{
	FILE *fp;
	char smess[SMAX];
	int l;
	if((fp=fopen(filename,"a"))!= NULL) {
		fprintf(fp,"\n%ld",dataset);
		fprintf(fp,"\t%s","All");
		fprintf(fp,"\t%.3f",r[nloc].totsites);
		fprintf(fp,"\t%.3f",r[nloc+1].totsites);
		fprintf(fp,"\t%.3f",r[nloc].bialsites);
		fprintf(fp,"\t%.3f",r[nloc+1].bialsites);
		fprintf(fp,"\t%.3f",r[nloc].multisites);
		fprintf(fp,"\t%.3f",r[nloc+1].multisites);
		fprintf(fp,"\t%.3f",r[nloc].sf);
		fprintf(fp,"\t%.3f",r[nloc+1].sf);
//		fprintf(fp,"\t%.3f",r[nloc].sfB);
//		fprintf(fp,"\t%.3f",r[nloc+1].sfB);
		fprintf(fp,"\t%.3f",r[nloc].sfout);
		fprintf(fp,"\t%.3f",r[nloc+1].sfout);
		fprintf(fp,"\t%.3f",r[nloc].sxA);
		fprintf(fp,"\t%.3f",r[nloc+1].sxA);
		fprintf(fp,"\t%.3f",r[nloc].sxB);
		fprintf(fp,"\t%.3f",r[nloc+1].sxB);
//		fprintf(fp,"\t%.3f",r[nloc].sxAfB);
//		fprintf(fp,"\t%.3f",r[nloc+1].sxAfB);
//		fprintf(fp,"\t%.3f",r[nloc].sxBfA);
//		fprintf(fp,"\t%.3f",r[nloc+1].sxBfA);
		fprintf(fp,"\t%.3f",r[nloc].ss);
		fprintf(fp,"\t%.3f",r[nloc+1].ss);
		fprintf(fp,"\t%.3f",r[nloc].Wald);
		fprintf(fp,"\t%.3f",r[nloc+1].Wald);
		fprintf(fp,"\t%.7f",r[nloc].piA);
		fprintf(fp,"\t%.7f",r[nloc].piB);
		fprintf(fp,"\t%.7f",r[nloc].thetaA);
		fprintf(fp,"\t%.7f",r[nloc].thetaB);
		fprintf(fp,"\t%.7f",r[nloc].DA);
		fprintf(fp,"\t%.7f",r[nloc].DB);
		fprintf(fp,"\t%.7f",r[nloc].dAB);
		fprintf(fp,"\t%.7f",r[nloc].dnAB);
		fprintf(fp,"\t%.7f",r[nloc].FST);
		for(l=0;l<nloc;l++) {	/*loop over loci*/
			fprintf(fp,"\t%d",l);
			fprintf(fp,"\t%d",(int) r[l].totsites);
			fprintf(fp,"\t%d",(int) r[l].bialsites);
			fprintf(fp,"\t%d",(int) r[l].multisites);
			fprintf(fp,"\t%d",(int) r[l].sf);
//			fprintf(fp,"\t%d",(int) r[l].sfB);
			fprintf(fp,"\t%d",(int) r[l].sfout);
			fprintf(fp,"\t%d",(int) r[l].sxA);
			fprintf(fp,"\t%d",(int) r[l].sxB);
//			fprintf(fp,"\t%d",(int) r[l].sxAfB);
//			fprintf(fp,"\t%d",(int) r[l].sxBfA);
			fprintf(fp,"\t%d",(int) r[l].ss);
			fprintf(fp,"\t%d",(int) r[l].Wald);
			fprintf(fp,"\t%.7f",r[l].piA);
			fprintf(fp,"\t%.7f",r[l].piB);
			fprintf(fp,"\t%.7f",r[l].thetaA);
			fprintf(fp,"\t%.7f",r[l].thetaB);
			fprintf(fp,"\t%.7f",r[l].DA);
			fprintf(fp,"\t%.7f",r[l].DB);
			fprintf(fp,"\t%.7f",r[l].dAB);
			fprintf(fp,"\t%.7f",r[l].dnAB);
			fprintf(fp,"\t%.7f",r[l].FST);
		}
		fclose(fp);
	}     /*if can open filename*/
	else {    /*if cannot open filename*/
		sprintf(smess,"\nError in  write_polyl :\n\tCannot open file %s",filename);
		write(ERRORFILE,smess);
		exit(1);
	}
}	/*end of write_polyl*/

/*****************************************************************************/

void initialize_write_polyl(char *filename,int nloc)
/*initialize the output file where writing statistics of shared and fixed polymorphism to outputfile one line for each replicate dataset*/
{
	FILE *fp;
	char smess[SMAX];
	int l;
	if((fp=fopen(filename,"a"))!= NULL) {
		fprintf(fp,"\n%s","dataset");
		fprintf(fp,"\t%s","loc");
		fprintf(fp,"\t%s","totsites_avg");
		fprintf(fp,"\t%s","totsites_std");
		fprintf(fp,"\t%s","bialsites_avg");
		fprintf(fp,"\t%s","bialsites_std");
		fprintf(fp,"\t%s","multisites_avg");
		fprintf(fp,"\t%s","multisites_std");
		fprintf(fp,"\t%s","sf_avg");
		fprintf(fp,"\t%s","sf_std");
//		fprintf(fp,"\t%s","sfB_avg");
//		fprintf(fp,"\t%s","sfB_std");
		fprintf(fp,"\t%s","sfout_avg");
		fprintf(fp,"\t%s","sfout_std");
		fprintf(fp,"\t%s","sxA_avg");
		fprintf(fp,"\t%s","sxA_std");
		fprintf(fp,"\t%s","sxB_avg");
		fprintf(fp,"\t%s","sxB_std");
//		fprintf(fp,"\t%s","sxAfB_avg");
//		fprintf(fp,"\t%s","sxAfB_std");
//		fprintf(fp,"\t%s","sxBfA_avg");
//		fprintf(fp,"\t%s","sxBfA_std");
		fprintf(fp,"\t%s","ss_avg");
		fprintf(fp,"\t%s","ss_std");
		fprintf(fp,"\t%s","Wald_avg");
		fprintf(fp,"\t%s","Wald_std");
		fprintf(fp,"\t%s","piA");
		fprintf(fp,"\t%s","piB");
		fprintf(fp,"\t%s","thetaA");
		fprintf(fp,"\t%s","thetaB");
		fprintf(fp,"\t%s","DtajA");
		fprintf(fp,"\t%s","DtajB");
		fprintf(fp,"\t%s","divAB");
		fprintf(fp,"\t%s","netdivAB");
		fprintf(fp,"\t%s","FST");
		for(l=0;l<nloc;l++) {	/*loop over loci*/
			fprintf(fp,"\t%s_%d","loc",l);
			fprintf(fp,"\t%s_%d","totsites",l);
			fprintf(fp,"\t%s_%d","bialsites",l);
			fprintf(fp,"\t%s_%d","multisites",l);
			fprintf(fp,"\t%s_%d","sf",l);
//			fprintf(fp,"\t%s_%d","sfB",l);
			fprintf(fp,"\t%s_%d","sfout",l);
			fprintf(fp,"\t%s_%d","sxA",l);
			fprintf(fp,"\t%s_%d","sxB",l);
//			fprintf(fp,"\t%s_%d","sxAfB",l);
//			fprintf(fp,"\t%s_%d","sxBfA",l);
			fprintf(fp,"\t%s_%d","ss",l);
			fprintf(fp,"\t%s_%d","Wald",l);
			fprintf(fp,"\t%s_%d","piA",l);
			fprintf(fp,"\t%s_%d","piB",l);
			fprintf(fp,"\t%s_%d","thetaA",l);
			fprintf(fp,"\t%s_%d","thetaB",l);
			fprintf(fp,"\t%s_%d","DtajA",l);
			fprintf(fp,"\t%s_%d","DtajB",l);
			fprintf(fp,"\t%s_%d","divAB",l);
			fprintf(fp,"\t%s_%d","netdivAB",l);
			fprintf(fp,"\t%s_%d","FST",l);
		}
		fclose(fp);
	}     /*if can open filename*/
	else {    /*if cannot open filename*/
		sprintf(smess,"\nError in  write_polyl :\n\tCannot open file %s",filename);
		write(ERRORFILE,smess);
		exit(1);
	}
}	/*end of write_polyl*/

/*****************************************************************************/

void write_ABCstat(char *filename,long dataset,int nloc,struct result_poly *r)
/*write statistics of shared and fixed polymorphism to outputfile one line for each replicate dataset*/
{
	FILE *fp;
	char smess[SMAX];
	int l;
	float corr_pi;

	if((fp=fopen(filename,"a"))!= NULL) {
		fprintf(fp,"\n%ld",dataset);
		/*fprintf(fp,"\t%.3f",r[nloc].totsites);
		fprintf(fp,"\t%.3f",r[nloc+1].totsites);*/
		fprintf(fp,"\t%.3f",r[nloc].bialsites);
		fprintf(fp,"\t%.3f",r[nloc+1].bialsites);
		fprintf(fp,"\t%.3f",r[nloc].sf);
		fprintf(fp,"\t%.3f",r[nloc+1].sf);
//		fprintf(fp,"\t%.3f",r[nloc].sfB);
//		fprintf(fp,"\t%.3f",r[nloc+1].sfB);
		fprintf(fp,"\t%.3f",r[nloc].sxA);
		fprintf(fp,"\t%.3f",r[nloc+1].sxA);
		fprintf(fp,"\t%.3f",r[nloc].sxB);
		fprintf(fp,"\t%.3f",r[nloc+1].sxB);
//		fprintf(fp,"\t%.3f",r[nloc].sxAfB);
//		fprintf(fp,"\t%.3f",r[nloc+1].sxAfB);
//		fprintf(fp,"\t%.3f",r[nloc].sxBfA);
//		fprintf(fp,"\t%.3f",r[nloc+1].sxBfA);
		fprintf(fp,"\t%.3f",r[nloc].ss);
		fprintf(fp,"\t%.3f",r[nloc+1].ss);
		fprintf(fp,"\t%.3f",r[nloc].Wald);
		fprintf(fp,"\t%.3f",r[nloc+1].Wald);
		fprintf(fp,"\t%.7f",r[nloc].piA);
		fprintf(fp,"\t%.7f",r[nloc+1].piA);
		fprintf(fp,"\t%.7f",r[nloc].piB);
		fprintf(fp,"\t%.7f",r[nloc+1].piB);
		fprintf(fp,"\t%.7f",r[nloc].thetaA);
		fprintf(fp,"\t%.7f",r[nloc+1].thetaA);
		fprintf(fp,"\t%.7f",r[nloc].thetaB);
		fprintf(fp,"\t%.7f",r[nloc+1].thetaB);
		fprintf(fp,"\t%.7f",r[nloc].DA);
		fprintf(fp,"\t%.7f",r[nloc+1].DA);
		fprintf(fp,"\t%.7f",r[nloc].DB);
		fprintf(fp,"\t%.7f",r[nloc+1].DB);
		fprintf(fp,"\t%.7f",r[nloc].dAB);
		fprintf(fp,"\t%.7f",r[nloc+1].dAB);
		fprintf(fp,"\t%.7f",r[nloc].dnAB);
		fprintf(fp,"\t%.7f",r[nloc+1].dnAB);
		fprintf(fp,"\t%.7f",r[nloc].FST);
		fprintf(fp,"\t%.7f",r[nloc+1].FST);
		corr_pi=pearson_corr_pi(r,nloc);
		fprintf(fp,"\t%.7f",corr_pi);
		fclose(fp);
	}     /*if can open filename*/
	else {    /*if cannot open filename*/
		sprintf(smess,"\nError in  write_ABCstat :\n\tCannot open file %s",filename);
		write(ERRORFILE,smess);
		exit(1);
	}
}	/*end of write_ABCstat*/

/*****************************************************************************/

void initialize_write_ABCstat(char *filename)
/*initialize the output file where writing statistics of shared and fixed polymorphism to outputfile one line for each replicate dataset*/
{
	FILE *fp;
	char smess[SMAX];
	int l;
	if((fp=fopen(filename,"w"))!= NULL) {
		fprintf(fp,"\n%s","dataset");
		/*fprintf(fp,"\t%s","totsites_avg");
		fprintf(fp,"\t%s","totsites_std");*/
		fprintf(fp,"\t%s","bialsites_avg");
		fprintf(fp,"\t%s","bialsites_std");
		fprintf(fp,"\t%s","sf_avg");
		fprintf(fp,"\t%s","sf_std");
//		fprintf(fp,"\t%s","sfB_avg");
//		fprintf(fp,"\t%s","sfB_std");
		fprintf(fp,"\t%s","sxA_avg");
		fprintf(fp,"\t%s","sxA_std");
		fprintf(fp,"\t%s","sxB_avg");
		fprintf(fp,"\t%s","sxB_std");
//		fprintf(fp,"\t%s","sxAfB_avg");
//		fprintf(fp,"\t%s","sxAfB_std");
//		fprintf(fp,"\t%s","sxBfA_avg");
//		fprintf(fp,"\t%s","sxBfA_std");
		fprintf(fp,"\t%s","ss_avg");
		fprintf(fp,"\t%s","ss_std");
		fprintf(fp,"\t%s","Wald_avg");
		fprintf(fp,"\t%s","Wald_std");
		fprintf(fp,"\t%s","piA_avg");
		fprintf(fp,"\t%s","piA_std");
		fprintf(fp,"\t%s","piB_avg");
		fprintf(fp,"\t%s","piB_std");
		fprintf(fp,"\t%s","thetaA_avg");
		fprintf(fp,"\t%s","thetaA_std");
		fprintf(fp,"\t%s","thetaB_avg");
		fprintf(fp,"\t%s","thetaB_std");
		fprintf(fp,"\t%s","DtajA_avg");
		fprintf(fp,"\t%s","DtajA_std");
		fprintf(fp,"\t%s","DtajB_avg");
		fprintf(fp,"\t%s","DtajB_std");
		fprintf(fp,"\t%s","divAB_avg");
		fprintf(fp,"\t%s","divAB_std");
		fprintf(fp,"\t%s","netdivAB_avg");
		fprintf(fp,"\t%s","netdivAB_std");
		fprintf(fp,"\t%s","FST_avg");
		fprintf(fp,"\t%s","FST_std");
		fprintf(fp,"\t%s","corr_pi");
		fclose(fp);
	}     /*if can open filename*/
	else {    /*if cannot open filename*/
		sprintf(smess,"\nError in  initialize_write_ABCstat :\n\tCannot open file %s",filename);
		write(ERRORFILE,smess);
		exit(1);
	}
}	/*end of initialize_write_ABCstat*/

/*****************************************************************************/
void write_stat_alleles(const char *filename,long gen,long na,float ne,float he)
/*writes value of naobs neobs heobs*/
{
	FILE *fp;
	char smess[SMAX];

	if((fp=fopen(filename,"a"))!= NULL) {
		if(gen == 1) fprintf(fp,"\nValues of statistics for each generation under equilibrium phase\n gen     \tnaobs\tneobs\theobs");
		fprintf(fp,"\n  %8ld\t%7ld\t%7.3f\t%7.3f",gen,na,ne,he);
		fclose(fp);
	}     /*if can open filename*/
	else {    /*if cannot open filename*/
		sprintf(smess,"\nError in write_stat_alleles :\n\tCannot open file %s",filename);
		write(ERRORFILE,smess);
		exit(1);
	}
}	  /*end of write_stat_alleles*/


/**********************************************************************************/

void write_allelic_freq(const char *filename,long gen,long na,struct typeallele *allele,
					float *allele_order, int model_type)
{
	FILE *fp;
	char smess[SMAX];
	long freqall,i,j;
	float id;

	if((fp=fopen(filename,"a"))== NULL) {
 		sprintf(smess,"\nError in write_stat_alleles :\n\tCannot open file %s",filename);
		write(ERRORFILE,smess);
		exit(1);
	}     /*if cannot open filename*/
	for(i=0;i<na;i++) {	  /*loop over alleles in allele_order*/
		if(model_type == DOMINANT) {
			id=allele_order[i];			
			for(j=0;j<na;j++) {	/*loop over alleles in allele to find id
											in order to record its order of dominance*/
				if(id==allele[j].id) {	/*found the right allele*/
					freqall=allele[j].freq;	/*gets frequency of allele i*/
					fprintf(fp,"\t%8ld",freqall);
					break;
				}	/*end of if found the right allele*/
			}	/*end loop  allele*/
		}	/*end of if(model_type == DOMINANT)*/
	}	/*end of for loop over alleles in allele_order*/
	fclose(fp);
}	  /*end of write_allelic_freq*/

/**********************************************************************************/

void writeTi(char *outputfile,int na,int nsimul,long **Ti)
{
	FILE *fp;
	int inodes;
	float mean,stdev;
	long count;

	if((fp=fopen(outputfile,"at"))==NULL) {
	   /*if cannot open filename*/
		write(ERRORFILE,"\ncannot open outputfile in writeTi");
		exit(1);
	} 
	fprintf(fp,"\n\tTi");
	for(inodes=na;inodes>1;inodes--) fprintf(fp,"\t");
	fprintf(fp,"Tistd");
	for(inodes=na;inodes>1;inodes--) fprintf(fp,"\t");
	fprintf(fp,"\nTree"); 
	for(inodes=na;inodes>1;inodes--) fprintf(fp,"\t%d",inodes);
	for(inodes=na;inodes>1;inodes--) fprintf(fp,"\t%d",inodes);
	/*for(itrees=0;itrees<nsimul;itrees++) {
		fprintf(fp,"\n%d",itrees+1);
		for(inodes=na;inodes>1;inodes--) fprintf(fp,"\t%7ld",Ti[inodes][itrees]); 
		for(inodes=na;inodes>1;inodes--) fprintf(fp,"\t%7ld",Ti[inodes][itrees]*inodes*(inodes-1)); 
    } */
	fprintf(fp,"\n%s","mean");
	for(inodes=na;inodes>1;inodes--) {
/*		avestdev_longint(Ti[inodes],nsimul,&mean,&stdev);*/
		avestd_long_miss(Ti[inodes],nsimul,MISSING,&mean,&stdev,&count);
		fprintf(fp,"\t%8.2f",mean);  
   	}
	for(inodes=na;inodes>1;inodes--) {
/*		avestdev_longint(Ti[inodes],nsimul,&mean,&stdev); */
		avestd_long_miss(Ti[inodes],nsimul,MISSING,&mean,&stdev,&count);
		fprintf(fp,"\t%8.2f",mean*inodes*(inodes-1));  
   	}
	fprintf(fp,"\n%s","SD");
	for(inodes=na;inodes>1;inodes--) {
/*		avestdev_longint(Ti[inodes],nsimul,&mean,&stdev); */
		avestd_long_miss(Ti[inodes],nsimul,MISSING,&mean,&stdev,&count);
		fprintf(fp,"\t%8.2f",stdev);  
   	}
	for(inodes=na;inodes>1;inodes--) {
/*		avestdev_longint(Ti[inodes],nsimul,&mean,&stdev); */
		avestd_long_miss(Ti[inodes],nsimul,MISSING,&mean,&stdev,&count);
		fprintf(fp,"\t%8.2f",stdev*inodes*(inodes-1));  
   	}
	fprintf(fp,"\n%s","trees");
	for(inodes=na;inodes>1;inodes--) {
/*		avestdev_longint(Ti[inodes],nsimul,&mean,&stdev);*/
		avestd_long_miss(Ti[inodes],nsimul,MISSING,&mean,&stdev,&count);
		fprintf(fp,"\t%8ld",count);  
   	}

   	fclose(fp);   
}	/*end of writeTi*/

/**********************************************************************************/

void writeTi_order(char *outputfile,int na,int nsimul,long **Ti_order)
{
	FILE *fp;
	int inodes;
	float mean,stdev;
	long count;

	if((fp=fopen(outputfile,"at"))==NULL) {
	   /*if cannot open filename*/
		write(ERRORFILE,"\ncannot open outputfile in writeTi_order");
		exit(1);
	} 
	fprintf(fp,"\n%s","dom_level");
	for(inodes=na;inodes>1;inodes--) {
/*		avestdev_longint(Ti[inodes],nsimul,&mean,&stdev);*/
		avestd_long_miss(Ti_order[inodes],nsimul,MISSING,&mean,&stdev,&count);
		fprintf(fp,"\t%8.2f",mean/10000);  
   	}
   	fclose(fp);   
}	/*end of writeTi_order*/



/******************************************************************************/

void write_life_allele(struct type_all *all,long gen,int model_type)
{

	FILE *fp;
	long i,k,maxfreq,maxfreq_i;
	char smess[SMAX];
	static int count=0;
	float mean,stdev;

	if(! count++) {	  /*if first call of this procedure*/
		if((fp=fopen(LIFE_ALLELE_FILE,"wt"))==NULL) {
		   /*if cannot open filename*/
			sprintf(smess,"cannot open file %s in write_life_allele",LIFE_ALLELE_FILE);
			write(ERRORFILE,smess);
			exit(1);
		}
		/*fprintf(fp,"allele_id\tlifetime\tmean_freq\tmax_freq");
		if (model_type == DOMINANT) 
			fprintf(fp,"\torder_atmax\tinitial order\tfinal order\tmean order");*/
		fprintf(fp,"allele_id\tlifetime\tmax_freq");
		if (model_type == DOMINANT) fprintf(fp,"\tinitial order\tfinal order");
		fclose(fp);
	}
	if((fp=fopen(LIFE_ALLELE_FILE,"at"))==NULL) {
	   /*if cannot open filename*/
		sprintf(smess,"cannot open file %s in write_life_allele",LIFE_ALLELE_FILE);
		write(ERRORFILE,smess);
		exit(1);
	}
	/*computes the number of generations since the origin of the allele*/
	k=gen - (all->time[(all->zero)-1])-1;	/*minus 1 because extinct in this generation*/
	fprintf(fp,"\n%10.3f\t%ld",all->id,k);
	/*avestdev_unsignedshort(all->freqtime, k+1,&mean,&stdev);*//*computes mean of freqtime across generations*/
	/*fprintf(fp,"\t%7.2f",mean);*/
	/*searches maximum of freqtime across generations*/
		/*maxfreq=0;
		for(i=0;i<=k;i++) {
			if(all->freqtime[i] > maxfreq) {
				maxfreq=all->freqtime[i];
				maxfreq_i=i; 
			}
		}*/
	maxfreq=all->maxfreq;
	fprintf(fp,"\t%ld",maxfreq);
	if (model_type == DOMINANT) {
		/*fprintf(fp,"\t%7.4f",all->order[maxfreq_i]);*/
		/*fprintf(fp,"\t%7.4f\t%7.4f",all->order_initial,all->order[k]);*/
		fprintf(fp,"\t%7.4f\t%7.4f",all->order_initial,all->order_actual);
		/*avestdev_float(all->order, k+1,&mean,&stdev);  computes mean of order across generations*/
		/*fprintf(fp,"\t%7.2f",mean);*/
	}
   	fclose(fp);   
}	/*end of write_life_allele*/


/****************************************************************************/

void avestd_long_miss(long int data[], unsigned long int n, int missing, float *ave, float *std,long *count_obs)
/*given array data[1..n], returns its mean as ave and its standard deviation as std
value in array data equal to missing are not taken into account*/
{
	unsigned long j;
	double s,ep,var,avetemp;
	unsigned long count=0;
	
	for (avetemp=0.0,j=0;j<n;j++) if(data[j] != missing) {   /*skip over missing values*/
		avetemp += (double) data[j];
		count++;
	}
	if(count);
	if(count) avetemp = avetemp/(double)count;
	*ave = (float) avetemp;
	var=ep=0.0;
	for (j=0;j<n;j++) if(data[j] != missing) {
		s=data[j]-(avetemp);
		ep += s;
		var += s*s;
	}
	if(count>1) var =(var-ep*ep/(double) count)/(double)(count-1.0);
	else var = 0.0;
	*std = (float) sqrt(var);
	*count_obs = count;
}  /*end of procedure avestd*/


/*****************************************************************************/

void write_neobs_int(const char *filename,long gen,float ne1,float ne2, long int_equilib)
/*writes average value of neobs over one interval*/
{
	FILE *fp;
	char smess[SMAX];

	if((fp=fopen(filename,"a"))!= NULL) {
		if(gen == int_equilib) fprintf(fp,"\nAverage values of neobs over each interval\n gen     \tneobs1\tneobs2");
		fprintf(fp,"\n  %8ld\t%7.3f\t%7.3f",gen,ne1,ne2);
		fclose(fp);
	}     /*if can open filename*/
	else {    /*if cannot open filename*/
		sprintf(smess,"\nError in write_neobs_int :\n\tCannot open file %s",filename);
		write(ERRORFILE,smess);
		exit(1);
	}
}	  /*end of write_neobs_int*/

