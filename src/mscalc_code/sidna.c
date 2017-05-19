#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spmain.h"
#include "sprun.h"
#include "nr.h"

void compute_polyl(int nseqA,int nseqB,int nsl,int nspolyl,
				   char **seqAlhs,char **seqBlhs,char *seqOls,struct result_poly *r)
/*computes the statistics of fixed and shared sites for locus l
	as well as nucleotide diversity and Watterson theta for each species
	as well as rough and net divergence between both species*/
{
	int h,h1,h2;
	int s;
	char first_base,current_site,outgroup;
	struct biallelic_sites biA,biB;
	int is_multiallelic;
	int count_bialsites=0;
	int count_multisites=0;
	int count_sf=0;
//	int count_sfB=0;
	int count_sfout=0;
	int count_sxA=0;
	int count_sxB=0;
//	int count_sxAfB=0;
//	int count_sxBfA=0;
	int count_ss=0;
	int	count_pair;
	double sumpairdif;
	int count_dif;
	int count_segr;
	float a,DTajima;
	float piT;
	int count_pair_piT;
	int *Walds;		/*vector of runs of fixed and polymorphic sites according to Miguel Navascues*/


	Walds = alloc_int_vector(nsl);
	for(s=0;s<nspolyl;s++) {	/*starts loop over polymorphic sites at locus l*/
		first_base=seqAlhs[0][s];
		/*is_biallelic =0;*/
		is_multiallelic=0;
		biA.n=1;
		biA.b1=first_base;
		for(h=1;h<nseqA;h++) {	/*starts loop over haplotypes at species A*/
			current_site=seqAlhs[h][s];
			if(current_site != first_base) {	/* if base different than first haplotype*/
				if(biA.n == 1) {				/* the site was not known as biallelic within A*/
					biA.n=2;
					biA.b2=current_site;
				} else if(biA.n == 2) {		/* the site was known as biallelic*/
					if(current_site != biA.b2) {	/*more than 2 different bases then multiallelic*/
						biA.n=3;
						is_multiallelic=1;
						break;
					}
				}	/*end else if*/
			}	/*end if current_site != first_base*/
		}	/*end loop over haplotypes at species A*/
		if(is_multiallelic == 1) {	/*multiple hits*/
			count_multisites++;
			continue;	/*skip to next site*/
		}
		first_base=seqBlhs[0][s];
		biB.n=1;
		biB.b1=first_base;
		for(h=1;h<nseqB;h++) {	/*starts loop over haplotypes at species B*/
			current_site=seqBlhs[h][s];
			if(current_site != first_base) {	/* if base different than first haplotype*/
				if(biB.n == 1) {				/* the site was not known as biallelic*/
					biB.n=2;
					biB.b2=current_site;
					/*is_biallelic=1;*/
				} else if(biB.n == 2) {		/* the site was known as biallelic*/
					if(current_site != biB.b2) {	/*more than 2 different bases then multiallelic*/
						biB.n=3;
						/*is_biallelic=0;*/
						is_multiallelic=1;
						break;
					}
				}	/*end else if*/
			}	/*end if current_site != first_base*/
		}	/*end loop over haplotypes at species B*/
		if(is_multiallelic == 1) {	/*multiple hits*/
			count_multisites++;
			continue;	/*skip to next site*/
		}
		outgroup=seqOls[s];
		/*now tests whether multiple hits and if not determines the type of polymorphism*/
		if(biA.n == 1) {	/*fixed in A*/
			if(biB.n == 1) {	/*fixed in B*/
				if(biA.b1 == biB.b1) {	/*same base fixed in A and B*/
					if(biA.b1 == outgroup)	continue; /*non polymorphic site skip to next site*/
					else {	/*site fixed along outgroup*/
						count_sfout++;
						count_bialsites++;
						continue;
					}	/*end of else*/
				} else {	/* different base fixed in A and B*/
					if(biA.b1 == outgroup) {	/* base in B is derived*/
						count_sf++;
						Walds[s]=1;
						count_bialsites++;
						continue;
					} else if(biB.b1 == outgroup) {	/* base in A is derived*/
						count_sf++;
						Walds[s]=1;
						count_bialsites++;
						continue;
					} else {	/* outgroup different than A and than B => multiallelic thus skip to next site*/
						count_multisites++;
						continue;	
					}
				}	/*end of else*/
			} else {	/*fixed in A but biallelic in B*/
				if((biA.b1 == biB.b1) || (biA.b1 == biB.b2)) {	/* fixed in A for one of the two bases in B*/
					if(outgroup == biA.b1) {	/*unique derived polymorphism in B*/
						count_sxB++;
						Walds[s]=0;
						count_bialsites++;
						continue;
					} else if((outgroup == biB.b1) || (outgroup == biB.b2)) {	/*unique ancestral polymorphism in B*/
						count_sxB++;
						Walds[s]=0;
						count_bialsites++;
						continue;
					} else {	/*outgroup different than 2 bases in B ==> multiallelic then next site*/
						count_multisites++;
						continue;	
					}
				} else {	/*base fixed in A is different than both in B*/
					count_multisites++;
					continue;	
				}
			}	/*end of else*/
		} else {	/*biallelic in A*/
			if(biB.n == 1) {	/*biallelic in A but fixed in B*/
				if((biB.b1 == biA.b1) || (biB.b1 == biA.b2)) {	/* fixed in B for one of the two bases in A*/
					if(outgroup == biB.b1) {	/*unique derived polymorphism in A*/
						count_sxA++;
						Walds[s]=0;
						count_bialsites++;
						continue;
					} else if((outgroup == biA.b1) || (outgroup == biA.b2)) {	/*unique ancestral polymorphism in A*/
						count_sxA++;
						Walds[s]=0;
						count_bialsites++;
						continue;
					} else {	/*outgroup different than 2 bases in A ==> multiallelic then skip to next site*/
						count_multisites++;
						continue;	
					}
				} else {	/*base fixed in B is different than both in A*/
					count_multisites++;
					continue;	
				}
			} else {	/*biallelic in both A and B*/
				if((biA.b1 == biB.b1) && (biA.b2 == biB.b2)) {	/*shared polymorphism*/
					if((outgroup == biA.b1) || (outgroup == biA.b2)) {	/*shared polym at biallelic site*/
						count_ss++;
						Walds[s]=0;
						count_bialsites++;
						continue;
					} else {	/*multiallelic*/
						count_multisites++;
						continue;	
					}
				} else if((biA.b1 == biB.b2) && (biA.b2 == biB.b1)) {	/*shared polymorphism*/
					if((outgroup == biA.b1) || (outgroup == biA.b2)) {	/*shared polym at biallelic site*/
						count_ss++;
						Walds[s]=0;
						count_bialsites++;
						continue;
					} else {	/*multiallelic because outgroup is different*/
						count_multisites++;
						continue;	
					}
				} else {	/* polymorphism for different pairs of bases in A and B*/
					count_multisites++;
					continue;	
				}
			}	/*end of else*/
		}	/*end of else*/
	}	/*end loop over sites*/
	r->totsites=nsl;
	r->bialsites=count_bialsites;
	r->multisites=count_multisites;
	r->sf=count_sf;
//	r->sfB=count_sfB;
	r->sfout=count_sfout;
	r->sxA=count_sxA;
	r->sxB=count_sxB;
//	r->sxAfB=count_sxAfB;
//	r->sxBfA=count_sxBfA;
	r->ss=count_ss;
	r->Wald=(float) WaldWolfowitz(nspolyl, Walds, 0);
	/*now computes pi for species A */
	count_pair=0;
	sumpairdif=0.0;
	for(h1=0;h1<(nseqA-1);h1++) {	/*starts loops over pairs of haplotypes at species A*/
		for(h2=h1+1;h2<nseqA;h2++) {
			count_dif=0;
			for(s=0;s<nspolyl;s++) {	/*starts loop over sites at locus l*/
				if(seqAlhs[h1][s] != seqAlhs[h2][s]) count_dif++;
			}
			count_pair++;
			sumpairdif += count_dif/(double)nsl;
				/*sumpairdif_sq += count * count;*/
		}
	}	/*end loops over pairs of haplotypes */
	r->piA=(float) sumpairdif/(float)count_pair;
	piT = sumpairdif;
	count_pair_piT = count_pair;
		/**avgpairdif_std= (double) (sumpairdif_sq - (countpair*(*avgpairdif)*(*avgpairdif)))/(countpair-1.0); 	*/
	/*now computes theta for species A */
	count_segr=0;
	for(s=0;s<nspolyl;s++) {	/*starts loop over sites at locus l*/
		first_base=seqAlhs[0][s];
		for(h=1;h<nseqA;h++) {	/*starts loop over haplotypes at species A*/
			if(seqAlhs[h][s] != first_base) {	/*current site is segregating*/
				count_segr++;
				break;
			}
		}	/*end loop over haplotypes in A*/
	}	/*end loop over sites*/
	a=0.0F;
	for(h=1;h<nseqA;h++) {	/*computes Watterson standardization parameter a*/
		a += (float) 1.0/(float) h; 
	}
	r->thetaA=(float) (count_segr/a)/(float) nsl;
	compute_DTajima(nseqA,count_segr,r->piA*nsl,&DTajima);	/*Warning Pi is multiplied by number of sites*/
	r->DA=DTajima;
	/*now computes pi for species B */
	count_pair=0;
	sumpairdif=0.0;
	for(h1=0;h1<(nseqB-1);h1++) {	/*starts loops over pairs of haplotypes at species B*/
		for(h2=h1+1;h2<nseqB;h2++) {
			count_dif=0;
			for(s=0;s<nspolyl;s++) {	/*starts loop over sites at locus l*/
				if(seqBlhs[h1][s] != seqBlhs[h2][s]) count_dif++;
			}
			count_pair++;
			sumpairdif += count_dif/(double)nsl;
				/*sumpairdif_sq += count * count;*/
		}
	}	/*end loops over pairs of haplotypes */
	r->piB=(float) sumpairdif/(float)count_pair;
	piT += sumpairdif;
	count_pair_piT += count_pair;
		/**avgpairdif_std= (double) (sumpairdif_sq - (countpair*(*avgpairdif)*(*avgpairdif)))/(countpair-1.0); 	*/
	/*now computes theta for species B */
	count_segr=0;
	for(s=0;s<nspolyl;s++) {	/*starts loop over sites at locus l*/
		first_base=seqBlhs[0][s];
		for(h=1;h<nseqB;h++) {	/*starts loop over haplotypes at species B*/
			if(seqBlhs[h][s] != first_base) {	/*current site is segregating*/
				count_segr++;
				break;
			}
		}	/*end loop over haplotypes in B*/
	}	/*end loop over sites*/
	a=0.0F;
	for(h=1;h<nseqB;h++) {	/*computes Watterson standardization parameter a*/
		a += (float) 1.0/(float) h; 
	}
	r->thetaB=(float) (count_segr/a)/(float) nsl;
	compute_DTajima(nseqB,count_segr,r->piB*nsl,&DTajima);
	r->DB=DTajima;
	/*now computes average pairwise divergence between species A and B*/
	count_pair=0;
	sumpairdif=0.0;
	for(h1=0;h1<nseqA;h1++) {	/*starts loop over haplotypes at species A*/
		for(h2=0;h2<nseqB;h2++) { /*starts loop over haplotypes at species B*/
			count_dif=0;
			for(s=0;s<nspolyl;s++) {	/*starts loop over sites at locus l*/
				if(seqAlhs[h1][s] != seqBlhs[h2][s]) count_dif++;
			}
			count_pair++;
			sumpairdif += count_dif/(double)nsl;
				/*sumpairdif_sq += count * count;*/
		}
	}	/*end loops over pairs of haplotypes */
	r->dAB=(float) sumpairdif/(float)count_pair;
	r->dnAB=r->dAB-(r->piA+r->piB)/2.0;
	piT += sumpairdif;
	count_pair_piT += count_pair;
	piT /= (float) count_pair_piT;
	if (piT < 1.0e-7) r->FST=MISSING;
	else r->FST=(piT-(r->piA+r->piB)/2.0)/piT;
	free(Walds);
}	/*end of compute_polyl*/

/**********************************************************************************/

/*From Miguel Navascues*/
// calculation of p-values has been removed

// Function to calculate number of runs or its standardize (i.e. normalized) value
// ARGUMENTS:
//   number_of_elements: length of the sequence
//   IOseq: pointer to a set of integers that represents the sequence to test
//          this sequence should contain exclusively "1" and "0" (ones and zeros)
//	 std: which statistic should be returned: 0=number of runs; 1=standardize number of runs
//        with any other value the number of runs is returned
//        (note that the number of runs is an integer but will be returned as a double)
double WaldWolfowitz(int number_of_elements, int *IOseq, int std)
{

	int i;
	int runs;
	int n0, n1, n;
	double mean;
	double variance;
	double z;

	if (std!=1 || std!=0) std=0;
	
	n=number_of_elements;

	//printf("Performing Wald-Wolfowitz test (aka Runs Test) on sequence\n");

	runs=1;
	n0=n1=0;
	for (i=0;i<number_of_elements;i++){
		if (IOseq[i]!=0 && IOseq[i]!=1){
			printf("sequence contains elements different from 1 and 0!\n");
			return -9;
		}
		if (IOseq[i]==0) n0++;
		if (IOseq[i]==1) n1++;
		if (i>0){
			if (IOseq[i]!=IOseq[i-1]) runs++;
		}
	}
	if (n0<5 || n1<5) return -9;
	if (n0+n1!=number_of_elements) printf("n0+n1<>n \n");
	if (runs>number_of_elements) printf("runs>n \n");

	if (std==1){
		mean = 1.0 + ( ((double)n0 * (double)n1 * 2.0) / (double)number_of_elements );
		variance = ( (mean-1) * (mean-2) ) / ((double)number_of_elements-1);
		z = ( (double)runs - mean ) / sqrt(variance);
		return z;
	}		
	if (std==0) return runs;
}


void compute_DTajima(int nseq,int nsegsites,float avgpairdif,float *D)
{
	float a1=0,a2=0,b1,b2,c1,c2,e1,e2;
	double Dtemp;
	int i;
		
	for(i=1;i<nseq;i++) {
		a1 += (float) 1.0/i; 
		a2 += (float) 1.0/i/i;
	}
	b1=(float) (nseq+1.0)/3.0/(nseq-1.0);
	b2=(float) 2.0*(nseq*nseq+nseq+3.0)/9.0/nseq/(nseq-1.0);
	c1=b1-(1.0/a1);
	c2=b2-((nseq+2.0)/a1/nseq)+(a2/a1/a1);
	e1=c1/a1;
	e2=c2/(a1*a1+a2);
	if(nsegsites == 0) {
		Dtemp=0.0F;	/*Warning put a value of 0 if no segregating sites!!!!!!*/
	} else	Dtemp=(double) (avgpairdif-((double) nsegsites/(double)a1))/sqrt(e1*nsegsites+e2*nsegsites*(nsegsites-1)); 
	*D=(float) Dtemp;
}	/*end of procedure compute_DTajima*/        

/**********************************************************************************************/



void create_anc_seq(int nnodes,int **site,int nsites,long int *seed) {
	int y,pickbase;
	for(y=0;y<nsites;y++) {		/*loop for different 
sites*/
		pickbase=4*ran1(seed);
		switch(pickbase) {
			case 0:
				site[nnodes-1] [y]='A';
				break;
			case 1:			        
				site[nnodes-1] [y]='G';
				break;
			case 2:			        
				site[nnodes-1] [y]='C';
				break;
			case 3:			        
				site[nnodes-1] [y]='T';
				break;
			default:
			printf("\nrandom generator problem in create_anc_seq"); 
			write_info("erreur","\n\nrandom generator problem in create_anc_seq"); exit(1);
		}	/*end of switch*/
	}	/*end of for sites loop*/
}	/*end procedure create_anc_seq*/
/**********************************************************************************/
void apply_mutation_toseq(int nmut,int *seq,int nsites, long int *seed)
{
	int i,picksite,basetomutate,pick;
	for(i=0;i<nmut;i++) {
		picksite=nsites*ran1(seed);
		basetomutate=seq[picksite];
		pick=3*ran1(seed);
		switch (basetomutate) {
			case 'A':
				switch (pick) { 
					case 0:
						seq[picksite]='G';
						break;
					case 1:
						seq[picksite]='C';
						break;
					case 2:
						seq[picksite]='T';
						break;
					default:
					printf("\nrandom generator problem in apply_mutation_toseq"); 
					write_info("erreur","\n\nrandom generator problem in apply_mutation_toseq"); exit(1);
				}	/*end of case where basetomutate=A*/
				break;
			case 'G':
				switch (pick) {
					case 0:
						seq[picksite]='A';
						break;
					case 1:
						seq[picksite]='C';
						break;
					case 2:
						seq[picksite]='T';
						break;
					default:
					printf("\nrandom generator problem in apply_mutation_toseq"); 
					write_info("erreur","\n\nrandom generator problem in apply_mutation_toseq"); exit(1);
				}	/*end of case where basetomutate=G*/
				break;
			case 'C':
				switch (pick) {
					case 0:
						seq[picksite]='A';
						break;
					case 1:
						seq[picksite]='G';
						break;
					case 2:
						seq[picksite]='T';
						break;
					default:
					printf("\nrandom generator problem in apply_mutation_toseq"); 
					write_info("erreur","\n\nrandom generator problem in apply_mutation_toseq"); exit(1);
				}	/*end of case where basetomutate=C*/
				break;
			case 'T':
				switch (pick) {
					case 0:
						seq[picksite]='A';
						break;
					case 1:
						seq[picksite]='G';
						break;
					case 2:
						seq[picksite]='C';
						break;
					default:
					printf("\nrandom generator problem in apply_mutation_toseq"); 
					write_info("erreur","\n\nrandom generator problem in apply_mutation_toseq"); exit(1);
				}	/*end of case where basetomutate=A*/
				break;
			default:
			printf("\nrandom generator problem in apply_mutation_toseq"); 
			write_info("erreur","\n\nrandom generator problem in apply_mutation_toseq");
			exit(1);
		}	/*end of switch basetomutate*/
	}	/*end of for different mutations to apply*/
}	/*end of procedure apply_mutation_toseq*/
/**********************************************************************************/

void build_sequences(struct typenode **tree,int nnodes,int nsample,int **site,int nsites,double mutrate2,long int *seed)    {
	int nodeid,i,y,ndesc,descid,picknmut;
	long int timeup;
	
	create_anc_seq(nnodes,site,nsites,seed); 
	for(nodeid=nnodes-1;nodeid>(nsample-1);nodeid--) {      /*loop ancestor nodes*/ 
		ndesc=tree[nodeid]->nnext;
		for(i=0;i<ndesc;i++) {          /*loop descendants of each ancestor nodes*/ 
			descid=tree[nodeid]->descendant[i]->mut_up;
			for(y=0;y<nsites;y++) site[descid] [y]=site[nodeid] [y]; 
			timeup=tree[descid]->time - tree[nodeid]->time; 
			picknmut=poidev(timeup*mutrate2,seed); 
			apply_mutation_toseq(picknmut,site[descid],nsites,seed);
		}	/*end loop descendants*/
	}	/*end loop ancestor nodes*/
}	/*end procedure build_sequences*/ 

/**********************************************************************************/
void count_segr_sites(int **site,int nsites,int nseq,int *nsegsites)
{
	int x,y,base,segr,count=0;
	for(y=0;y<nsites;y++) {         /* loop different sites*/
		base=site[0][y];
		segr=0;
		for(x=1;x<nseq;x++) {   /*loop sequences*/
			if(site[x][y] != base) {        /*if segregating site*/
				segr=1;
				break;
			}	        /*if segregating site*/
		}	/*loop sequences*/
		if(segr) count++;
	}	/*end loop different sites*/
	*nsegsites=count;
}	/*end procedure count_segr_sites*/ 

/**********************************************************************************/
void comp_avgpair_nucldiff(int **site,int nsites,int nseq,double *avgpairdif,double *avgpairdif_std) {
	long int sumpairdif=0,sumpairdif_sq=0,countpair=0;
	int i,j,y,count;
	for(i=0;i<(nseq-1);i++) {
		for(j=1;j<nseq;j++) {
			count=0;
			for(y=0;y<nsites;y++) {
				if(site[i][y] != site[j][y]) count++;
			}
			countpair++;
			sumpairdif += count;
			sumpairdif_sq += count * count;
		}	/*end loop j*/
	}	/*end loop i*/
	*avgpairdif= (double) sumpairdif/countpair;
	*avgpairdif_std= (double) (sumpairdif_sq - (countpair*(*avgpairdif)*(*avgpairdif)))/(countpair-1.0); }	/*end procedure 
comp_avgpair_nucldiff*/
/**********************************************************************************/
void Tajima_test(int **site,int nsites,int nseq,float *D,int *nsgsites,float *avgpd,float *avgpd_std)
{
	float a1=0,a2=0,b1,b2,c1,c2,e1,e2;
	double avgpairdif,avgpairdif_std,Dtemp;
	int i,nsegsites;
		
	count_segr_sites(site,nsites,nseq,&nsegsites); 
	comp_avgpair_nucldiff(site,nsites,nseq,&avgpairdif,&avgpairdif_std);
	for(i=1;i<nseq;i++) {
		a1 += (float) 1.0/i; 
		a2 += (float) 1.0/i/i;
	}
	b1=(float) (nseq+1.0)/3.0/(nseq-1.0);
	b2=(float) 2.0*(nseq*nseq+nseq+3.0)/9.0/nseq/(nseq-1.0);
	c1=b1-(1.0/a1);
	c2=b2-((nseq+2.0)/a1/nseq)+(a2/a1/a1);
	e1=c1/a1;
	e2=c2/(a1*a1+a2);
	Dtemp=(double) (avgpairdif-((double) nsegsites/a1))/sqrt(e1*nsegsites+e2*nsegsites*(nsegsites-1)); 
	*D=Dtemp;
	*nsgsites=nsegsites;
	*avgpd=avgpairdif;
	*avgpd_std=avgpairdif_std;
}	/*end of procedure Tajima_test*/        
/**********************************************************************************/

float pearson_corr_pi(struct result_poly *resl,int nloc)
{
    int j;
    float yt,xt;
    float syy=0.0, sxy=0.0, sxx=0.0, ay=0.0, ax=0.0;

	char smess[SMAX];
    for (j = 0 ; j < nloc ; j++)
    {
/*sprintf(smess,"\n%.6f\t%.6f",resl[j].piA,resl[j].piB);
write(ERRORFILE,smess);*/
	}

    for (j = 0 ; j < nloc ; j++)
    {
        ax += resl[j].piA;
        ay += resl[j].piB;
    }
    ax /= nloc;
    ay /= nloc;

    for (j = 0 ; j < nloc ; j++)
    {
        xt=resl[j].piA-ax;
        yt=resl[j].piB-ay;
        sxx += xt*xt;
        syy += yt*yt;
        sxy += xt*yt;
    }
    return sxy/sqrt(sxx*syy);
}




