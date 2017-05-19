
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "spmain.h"

#define NR_END 1
#define FREE_ARG char*
extern long glob_tos_node_stack;/*global variable with index of next open stack location*/
extern long glob_max_node_stack;/*global variable with maximum number of  nodes in the stack*/
extern struct type_all **glob_node_stack;	/*array of pointers to nodes in the stack*/
extern long glob_min_tos_node_stack;/*glob var with minimum value taken by tos*/


struct result_poly *alloc_result_poly(long size)
/*allocates a vector of elements of structure type result_poly*/
{
	struct result_poly *p;
	
	p=(struct result_poly *) malloc((size_t) size * sizeof(struct result_poly)); 
	if(!p) { 
		write(ERRORFILE,"\nOut of memory in alloc_result_poly"); 
		exit(1); 
	}
	return p;
}

/**************************************************************************************/ 



struct typeallele *alloc_allele(long size)
/*allocates an array of elements of the structure type typeallele 
	for recording info about allele ID's and frequency in current population*/
{
	struct typeallele *p;
	
	p=(struct typeallele *) malloc((size_t) size * sizeof(struct typeallele)); 
	if(!p) { 
		write(ERRORFILE,"\nOut of memory in alloc_allele"); 
		exit(1); 
	}
	return p;
}

/**************************************************************************************/ 
struct typeallele **alloc_matrix_allele(long nrow,long ncol)
/*allocates an array of pointer to elements of the structure type typeallele 
	for recording info about allele ID's and frequency in current population*/
{
	struct typeallele **m;
	long i;

	/* allocate pointers to rows */
	m=(struct typeallele **) malloc((size_t) nrow * sizeof(struct typeallele*)); 
	if(!m) { 
		write(ERRORFILE,"\nOut of memory in alloc_matrix_allele"); 
		exit(1); 
	}
	/* allocate rows and set pointers to them */
	m[0]=(struct typeallele *) malloc((size_t) ((nrow*ncol)* sizeof(struct typeallele))); 
	if (!m[0]) { 
		write(ERRORFILE,"\nout of memory in alloc_matrix_allele"); 
		exit(1); 
	}
	for(i=1;i<nrow;i++) m[i]=m[i-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}	/*end of alloc_matrix_allele*/

/**************************************************************************************/ 

struct type_all **alloc_all(long size)
/*allocates memory to an array of pointers to struct of type type_all
  attention this does not allocates memory to the element of the structure
  these are allocated to by procedure alloc_nodes_to_all*/
{
	struct type_all **p;
	p=(struct type_all **) malloc((size_t) size * sizeof(struct type_all *));
	if(!p) { 
		write(ERRORFILE,"\nOut of memory in alloc_all");
		exit(1); 
	}
	return p;
}

/******************************************************************/

void alloc_nodes_to_all(struct type_all **all,long size) 
{
	struct type_all *node;
	long i;
	for(i=0;i<size;++i) {                     
		/*create size nodes pointed to by all */
		node= (struct type_all *)malloc((size_t) sizeof(struct type_all));
		if(!node) { 
			write(ERRORFILE,"\nOut of memory in alloc_nodes_to_all"); 
			exit(1);
		}
		all[i]=node;
	}
}    /* end of alloc_nodes_to_all */ 

/************************************************************************************/

void push_node_stack(struct type_all *node)
/*put an element on the stack p555 from book C the complete reference, by Schildt*/
{

	if(glob_tos_node_stack >= glob_max_node_stack) {
		write(ERRORFILE,"\nnode_stack is full");
		exit(1);
	}
	glob_node_stack[glob_tos_node_stack] = node;   
	glob_tos_node_stack++;
}
/***************************************************************************************/ 

struct type_all *pop_node_stack(void)
/*retrieve the top element from the stack*/
{
	glob_tos_node_stack--;
	if(glob_tos_node_stack <0) {
		write(ERRORFILE,"\nError in pop_node_stack: node_stack underflow");
		write(ERRORFILE,"\n\tthis means that you have to increase the value of glob_max_node_stack in spmain.c");
		/*increase_node_stack(INCREASE_STACK);*/
		exit(1);
	}
	if(glob_tos_node_stack < glob_min_tos_node_stack) glob_min_tos_node_stack = glob_tos_node_stack;   
	return glob_node_stack[glob_tos_node_stack];  
}

 /******************************************************************/ 

struct type_allele_distrib *alloc_allele_distrib(long size)
/*allocates an array of elements of the structure type type_allele_distrib 
	for recording info about allele frequency distribution over all replicate simulations*/
{
	struct type_allele_distrib *p;
	
	p=(struct type_allele_distrib *) malloc((size_t) size * sizeof(struct type_allele_distrib)); 
	if(!p) { 
		write(ERRORFILE,"\nOut of memory in alloc_allele_distrib"); 
		exit(1); 
	}
	return p;
}
 /******************************************************************/ 

struct type_allele_distrib_class *alloc_allele_distrib_class(long size)
/*allocates an array of elements of the structure type type_allele_distrib_class 
	for recording info about allele frequency distribution over all replicate simulations*/
{
	struct type_allele_distrib_class *p;
	
	p=(struct type_allele_distrib_class *) malloc((size_t) size * sizeof(struct type_allele_distrib_class)); 
	if(!p) { 
		write(ERRORFILE,"\nOut of memory in alloc_allele_distrib_class"); 
		exit(1); 
	}
	return p;
}

/******************************************************************/

struct type_tree *alloc_tree(long size)
{
	struct type_tree *p;
	/*if(DEBUGMEM) printf("\nprocedure alloc_tree\n");*/
	p=(struct type_tree *) malloc((size_t) size * sizeof(struct type_tree));
	if(!p) { 
		write(ERRORFILE,"\nout of memory in alloc_tree");
		exit(1); 
	}
	return p;
}  

/******************************************************************/ 
/******************************************************************/ 
/*	now standard vector and array allocation procedures 
	(adapted from Numerical recipes in C)*/
/******************************************************************/ 
/******************************************************************/ 

int *alloc_int_vector(long size)
/*allocate an int vector of size "size" will be defined as an array from 0 to size-1 */
{
	int *v;
	v=(int*)malloc((size_t) (size*sizeof(int)));
	if(!v) { 
		write(ERRORFILE,"\nout of memory in alloc_int_vector");
		exit(1); }
	return v;
}
 
/******************************************************************/ 

long *alloc_longint_vector(long size)
/*allocate a long int vector of size "size" =>0 to size-1*/ 
{ 
	long *v;
	v=(long *)malloc((size_t) (size*sizeof(long )));
	if(!v) { 
		write(ERRORFILE,"\nout of memory in alloc_longint_vector");
		exit(1); }
	return v;
}

/******************************************************************/

int *ivector(long nl,long nh)
/*allocate an int vector with subscript range v[nl..nh]*/
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if(!v) {
		write(ERRORFILE,"\nout of memory in ivector");
		exit(1);
	}
	return v-nl+NR_END;
}

void free_ivector(int *v, long nl, long nh)
/*free an int vector allocated with ivector*/
{
	free((FREE_ARG) (v+nl-NR_END));
}
/******************************************************************/





long **longmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a long int matrix with subscript range m[nrl..nrh][ncl..nch]
from Numerical Recipes in C */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	long  **m;
	/* allocate pointers to rows */
	m=(long  **) malloc((size_t)((nrow+NR_END)*sizeof(long *)));
	if (!m)  { 
		write(ERRORFILE,"\nout of memory in longmatrix");
		exit(1); 
	}
	m += NR_END;
	m -= nrl;
	/* allocate rows and set pointers to them */
	m[nrl]=(long *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(long)));
	if (!m[nrl]) { 
		write(ERRORFILE,"\nout of memory in longmatrix");
		exit(1); 
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for(i=nrl+1;i<=nrh;i++) {
		m[i]=m[i-1]+ncol;
	}
	/* return pointer to array of pointers to rows */
	return m;
} 

/***************************************************************************************/
 
void free_longmatrix(long **m, long nrl, long nrh, long ncl, long nch)
/* free a long int matrix allocated by longmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}
/******************************************************************/
float **floatmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;
	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m)  { 
		write(ERRORFILE,"\nout of memory in floatmatrix");
		exit(1); 
	}
	m += NR_END;
	m -= nrl;
	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float))); 
	if (!m[nrl]) { 
		write(ERRORFILE,"\nout of memory in floatmatrix"); 
		exit(1); 
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

/***************************************************************************************/

void free_floatmatrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by floatmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

/******************************************************************/

int **intmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;
	/* allocate pointers to rows */
	m=( int **) malloc((size_t)((nrow+NR_END)*sizeof( int*)));
	if (!m)  { 
		write(ERRORFILE,"\nout of memory in intmatrix");
		exit(1); 
	}
	m += NR_END;
	m -= nrl;
	/* allocate rows and set pointers to them */
	m[nrl]=( int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof( int))); 
	if (!m[nrl]) { 
		write(ERRORFILE,"\nout of memory in intmatrix"); 
		exit(1); 
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
} /***************************************************************************************/ 
void free_intmatrix( int **m, long nrl, long nrh, long ncl, long nch)
/* free an  int matrix allocated by longmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
} 


/*********************************************************************/

long *lvector(long nl, long nh)
/*allocate a long vector with subscript range[nl..nh]*/
{
	long *v;
	v=(long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long))); 
	if (!v) {
		printf("\nallocation failure in lvector()");
		exit(1);
	}
	return v-nl+NR_END;
}

/*********************************************************************/

void free_lvector(long  *v, long nl, long nh)
/*free a long  vector allocated with vector()*/
{
	free((FREE_ARG) (v+nl-NR_END));
}

/********************************************************************************/

float *vector(long int nl, long int nh)
/*allocate a float vector with subscript range[nl..nh]*/
{
	float *v;
v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float))); if (!v) {
		write(ERRORFILE,"\nAllocation failure in vector()");
		exit(1);
	}
	return v-nl+NR_END;
}
/*******************************************************************
**/
void free_vector(float *v, long int nl, long int nh)
/*free a float vector allocated with vector()*/
{
	free((FREE_ARG) (v+nl-NR_END));
}
/*******************************************************************
**/
double *dvector(long int nl, long int nh)
/*allocate a double vector with subscript range[nl..nh]*/
{
	double *v;
v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double))); if 
(!v) {
		write(ERRORFILE,"\nallocation failure in dvector()");
		exit(1);
	}
	return v-nl+NR_END;
}
/*******************************************************************
**/
void free_dvector(double *v, long int nl, long int nh)
/*free a double vector allocated with vector()*/
{
	free((FREE_ARG) (v+nl-NR_END));
}
/*******************************************************************
**/

/****************************************************************************/

int ***i3tensor(long nrl, long nrh,long ncl, long nch,long ndl,long ndh)
/*allocate a int 3tensor with subscript range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	int ***t;
	
	/*allocate pointers to pointers to rows*/
	t=(int ***) malloc((size_t)((nrow+NR_END)*sizeof(int **)));
	if (!t) {
		write("erreur","\nallocation failure 1 in f3tensor()");
		exit(1);
	}
	t += NR_END;
	t -= nrl;
	
	/*allocate pointers to rows and set pointers to them*/
	t[nrl]=(int **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int *)));
	if (!t[nrl]) {
		write("erreur","\nallocation failure 2 in f3tensor()");
		exit(1);
	}
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	
	/*allocate rows and set pointers to them*/
	t[nrl][ncl]=(int *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(int )));
	if (!t[nrl][ncl]) {
		write("erreur","\nallocation failure 3 in f3tensor()");
		exit(1);
	}
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}
	
	/*return pointer to array of pointers to rows*/
	return t;
}

void free_i3tensor(int ***t,long nrl, long nrh,long ncl, long nch,long ndl,long ndh)
/*free an int f3tensor allocated by f3tensor()*/
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}	

/****************************************************************************/

/****************************************************************************/

char ***char3tensor(long nrl, long nrh,long ncl, long nch,long ndl,long ndh)
/*allocate a char 3tensor with subscript range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	char ***t;
	
	/*allocate pointers to pointers to rows*/
	t=(char ***) malloc((size_t)((nrow+NR_END)*sizeof(char **)));
	if (!t) {
		write("erreur","\nallocation failure 1 in char3tensor()");
		exit(1);
	}
	t += NR_END;
	t -= nrl;
	
	/*allocate pointers to rows and set pointers to them*/
	t[nrl]=(char **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(char *)));
	if (!t[nrl]) {
		write("erreur","\nallocation failure 2 in char3tensor()");
		exit(1);
	}
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	
	/*allocate rows and set pointers to them*/
	t[nrl][ncl]=(char *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(char )));
	if (!t[nrl][ncl]) {
		write("erreur","\nallocation failure 3 in char3tensor()");
		exit(1);
	}
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}
	
	/*return pointer to array of pointers to rows*/
	return t;
}

void free_char3tensor(int ***t,long nrl, long nrh,long ncl, long nch,long ndl,long ndh)
/*free an int char3tensor allocated by f3tensor()*/
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}	


/******************************************************************/

char **charmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate an char matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	char **m;
	/* allocate pointers to rows */
	m=( char **) malloc((size_t)((nrow+NR_END)*sizeof( char*)));
	if (!m)  { 
		write(ERRORFILE,"\nout of memory in charmatrix");
		exit(1); 
	}
	m += NR_END;
	m -= nrl;
	/* allocate rows and set pocharers to them */
	m[nrl]=( char *) malloc((size_t)((nrow*ncol+NR_END)*sizeof( char))); 
	if (!m[nrl]) { 
		write(ERRORFILE,"\nout of memory in charmatrix"); 
		exit(1); 
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
} /***************************************************************************************/ 
void free_charmatrix( char **m, long nrl, long nrh, long ncl, long nch)
/* free an  char matrix allocated by longmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
} 


