#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "spmain.h"

#define NR_END 1
#define FREE_ARG char*
#ifndef ERRORFILE
	#define ERRORFILE "error.txt"
#endif


/*************************************************************************/

double sqr(float x)
{
	return (double) x*(double)x;
}

/*************************************************************************/

int roundxa(float x)
{
	int i ;
	if((x - (int) x)>= 0.5) i= ((int) x )+1;
	else i=(int) x;
	return i;
}

/*************************************************************************/

long roundxa_long(float x)
{
	long i ;
	if((x - (long) x)>= 0.5) i= ((long) x )+1;
	else i=(long) x;
	return i;
}

/*************************************************************************/

void get_time_short(char *smess)
/*get the current time and put only the day number and time in smess*/
{
	char s1[50],*s3;
	time_t t_current;
	struct tm *ptr;
	
	t_current = time(NULL);
	ptr = localtime(&t_current);
	sprintf(s1,"%s",asctime(ptr));
	s3=&s1[0];
	s3 = s3+7;
	s3[13]='\0'; 
    strcpy(smess,s3);
}

/***************************************************************************/

void get_time_full(char *smess)
/*get the current time and put it in full in smess*/
{
	char s1[50];
	time_t t_current;
	struct tm *ptr;
	
	t_current = time(NULL);
	ptr = localtime(&t_current);
	sprintf(s1,"%s",asctime(ptr));
    strcpy(smess,s1);
}

/****************************************************************************/

void write_info(char *filename,char *info)
/*  writes string info to logfile*/
{
	FILE *fp;
	
	if((fp=fopen(filename,"a"))!= NULL) {
		fprintf(fp," %s",info);
		fclose(fp);
	}     /*if can open filename*/
	else {    /*if cannot open filename*/
		printf("\nError in writing %s : cannot open file",info); exit(1);
	}
}	  /*end of write_info*/

/****************************************************************************/

void write(char *filename,char *info)
/* print string info to display and writes string info to file*/
{
	FILE *fp;
	
	printf("%s",info);
	if((fp=fopen(filename,"a"))!= NULL) {
		fprintf(fp,"%s",info);  
		fclose(fp);
	}     /*if can open filename*/
	else {    /*if cannot open filename*/
		printf("\nError in writing %s : cannot open file",info); 
		exit(1);
	}
}	  /*end of write*/


/****************************************************************************/

void write_onlytofile(char *filename,char *info)
/* print string info to display and writes string info to file*/
{
	FILE *fp;
	
/*	printf("%s",info);*/
	if((fp=fopen(filename,"a"))!= NULL) {
		fprintf(fp,"%s",info);  
		fclose(fp);
	}     /*if can open filename*/
	else {    /*if cannot open filename*/
		printf("\nError in writing %s : cannot open file",info); 
		exit(1);
	}
}	  /*end of write_onlytofile*/


/********************************************************************************/

void avestdev_longint(long int data[], unsigned long int n,float *ave,float *stdev)
/*given array data[0..n-1], returns its mean as ave and standard deviation(n-1) as stdev from Numerical recipes*/ {
		unsigned long int j;
		double s,ep,avetemp,var;
		
		for(avetemp=0.0,j=0;j<n;j++) avetemp += (double) data[j]; 
		if(n) avetemp /= n;
		*ave = (float) avetemp;
		var=ep=0.0;
		for (j=0;j<n;j++) {
			s=data[j]-(avetemp);
			ep += s;
			var += s*s;
		}
		if(n && (n-1)) var=(var-ep*ep/(double)n)/(double)(n-1.0);
		*stdev=sqrt(var);
}		/*end of function avestdev*/

/**********************************************************************/

void avestdev_float(float data[], unsigned long int n,float *ave,float *stdev)
/*given array data[0..n-1], returns its mean as ave and standard deviation(n-1) as stdev from Numerical recipes*/ {
		unsigned long int j;
		float s,ep,var;
		for(*ave=0.0F,j=0;j<n;j++) *ave += data[j]; 
		if(n) *ave /= n;
		var=ep=0.0;
		for (j=0;j<n;j++) {
			s=data[j]-(*ave);
			ep += s;
			var += s*s;
		}
		if(n && (n-1)) var=(var-ep*ep/(float) n)/(float) (n-1.0F);
		*stdev=sqrt(var);
}		/*end of function avestdev*/ 


/****************************************************************************/

void avestd_float_missing(float data[], unsigned long n, float missing, float *ave, float *std)
/*given array data[0..n-1], returns its mean as ave and its standard deviation as std
value in array data equal to missing are not taken into account*/
{
	unsigned long j;
	float s,ep,var;
	unsigned long count=0;
	
	for (*ave=0.0F,j=0;j<n;j++) if(data[j] != missing) {   /*skip over missing values*/
		*ave += data[j];
		count++;
	}
	if(count) *ave /= (float) count;
	else *ave=0.0F;
	var=ep=0.0F;
	for (j=0;j<n;j++) if(data[j] != missing) {
		s=data[j]-(*ave);
		ep += s;
		var += s*s;
	}
	if(count>1) var =(var-ep*ep/(float) count)/(float) (count-1.0F);
	else var = 0.0F;
	*std = (float) sqrt(var);
}  /*end of procedure avestd_missing*/




/********************************************************************************/

void avestdev_unsignedshort(unsigned short data[], unsigned long int n,float *ave,float *stdev)
/*given array data[0..n-1], returns its mean as ave and standard deviation(n-1) as stdev from Numerical recipes*/ {
		unsigned long int j;
		double s,ep,avetemp,var;
		
		for(avetemp=0.0,j=0;j<n;j++) avetemp += (double) data[j]; 
		if(n) avetemp /= n;
		*ave = (float) avetemp;
		var=ep=0.0;
		for (j=0;j<n;j++) {
			s=data[j]-(avetemp);
			ep += s;
			var += s*s;
		}
		if(n && (n-1)) var=(var-ep*ep/(double)n)/(double)(n-1.0);
		*stdev=sqrt(var);
}		/*end of function avestdev*/


/**********************************************************************/

double comp_std(long n,double sy,double sy2)
	/*computes de standard deviation from sum and sum of square of y with N=n, 
		p.53,55 Sokal and Rholf*/ 
{
	double ssq;

	if(n<2) return 0.0;
	ssq=sy2- ((sy*sy)/(double) n);
	return sqrt(ssq/(double) (n-1.0));  
}	        /*end of procedure comp_std*/ 

/*********************************************************************/
