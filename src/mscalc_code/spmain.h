
#define MAX_NLOC 25000	/*max number of loci in simulation datasets*/
#define SMAX 1000			/*size of string arrays used for messages*/



/*do not change the four next definitions*/
#define NEUTRAL 1			  /*neutral locus*/
#define GAMETOPHYTIC 2		  /*gametophytic SI*/
#define SSICOD 3				/*SSI with codominance in pollen and style*/
#define SSIDOMCOD 4				/*SSI with dom in pollen and codom in style*/
#define SSIDOMCOD_ONESTEP 5		/*Idem SSIdomcod but mutations only up or down one step*/
#define SSIDOM 6				/*SSI with dom in pollen and style*/
#define SSIDOM_ONESTEP 7		/*Idem SSIdom but mutations only up or down one step*/
#define SSIDOMCOD_INDEP 8		/*Idem SSIdomcod but mutation to absolute level*/
#define SSIDOM_INDEP 9			/*idem SSIdom but mutation to absolute level*/
#define OVERDOMINANT_S 10	/*model of symmetric overdominance with ovd_s as selection coefficient*/
#define HAPLOID_ONE 11	/*model of reproduction haploids with two loci*/
#define HAPLOID_TWO 12	/*model of reproduction haploids with one locus*/


#define MAX_VALUE_MODELS 12	/*max value of model_of_reproduction*/ 

/*#define POLLEN_CODOMINANCE 3  codominance in pollen*/
/*#define POLLEN_HIERARCHY 4	  linear dominance hierarchy in pollen*/
/*#define POLLEN_HIERARCHY_INDEP 5  dominance hierarchy in pollen indenpendant of current alleles*/
/*#define POLSTIGMA_HIERARCHY 6	linear dominance hierarchy in both stigma and pollen*/
/*you should set the variable model_of_reproduction in inputfile to the desired model of reproduction*/

#define NORMAL 0		/*model_type is set to NORMAL if no dominance relationships among alleles*/
#define DOMINANT 1		/*model_type is set to DOMINANT if dominance in hierarchical dominance classes*/
#define INDEPENDENT 2	/*model_type is set to INDEPENDENT if dominance in absolute value*/


#define MAX_DESC  10			/*max number of descendants per individual in struct typenode*/
#define NMAX 1000000          /*max number of individuals in the population*/ 
#define NSIMULMAX 60000     /*max number of simulation runs*/
#define MAX_MUT 5000      /*max number of items in array all.mut and all.time
							Attention should be a integer (<30000)*/
#define ERRORFILE "error.txt"  /*file to use for error messages*/
#define LIFE_ALLELE_FILE "splife.txt"  /*file to use for recording statistics about individual alleles*/
#define MAXGENEQ 100000	/*maximum number of generations to perform in waiting for an equilibrium
						of allelic frequencies*/
#define MAXGENGENEALOGY	 9000000 /*maximum number of generations to perform before ending the genealogical phase*/
/*#define MAX_AGE_ALLELE	 10000*/ 	/*maximum age of an individual allele, used for recording frequency and order changes*/
#define	MISSING -9999			/*value to use as missing value indicator*/
#define NALL_LIFE_MAX 100000	/*max number of alleles that can be followed for computation of probab of invasion/exit*/
#define NCLAS_LIFE_MAX 1000		/*max number of dominance classes for computation of probab of invasion/exit*/


/*defines the structure used to store results from shared sites, fixed sites per locus*/
struct result_poly{
	float totsites;
	float bialsites;
	float multisites;
	float sf;
//	float sfB;
	float sfout;
	float sxA;
	float sxB;
//	float sxAfB;
//	float sxBfA;
	float ss;
	float piA;
	float piB;
	float thetaA;
	float thetaB;
	float DA;
	float DB;
	float dAB;
	float dnAB;
	float FST;
	float Wald;
};

/*defines the structure used to store results from shared sites, fixed sites per locus*/
struct biallelic_sites{
	int n;
	char b1;
	char b2;
};

/*defines the structure used when inventories of alleles in the current population are made*/
struct typeallele{
	float id;		 /*identification of the allele = allele ID*/
	long freq;       /*number of occurence of the same id in the population*/
};

/*defines the structure used to implement the genealogical information, 
  following the method by Takahata & Nei,1990*/ 
struct type_all {
	float id;		/*allele ID*/
	int zero;      	/*position of first zero: =1 if first zero at mut[1]*/
	int nmut;           /*number of times this allele has mutated to new alleles*/
	int mut[MAX_MUT];	/*array with mutation id, according to Takahata & Nei*/
	long time[MAX_MUT];	/*array with time of occurence of mutations*/
	float order_par[MAX_MUT]; /*array with relative order of dominance of the parental alleles at time of mutation*/
	float order_initial;/*relative order of dominance of the allele when originated*/
	float order_actual;	/*current relative order of dominance of the allele*/
	int maxfreq;	/*maximum frequency of the allele during its life span*/
	/*unsigned short freqtime[MAX_AGE_ALLELE];*/ /*array of frequency of the allele at each subsequent generations*/ 
	/*float order[MAX_AGE_ALLELE];*/	/*array of relative order of the allele at each subsequent generation*/
};

/*defines the structure used to record allelic frequencies distribution over replicate runs*/
struct type_allele_distrib {
	unsigned int count;		/*number of occurence of allele with given freq in all simulations*/
	float rel_order;		/*order in allele_order if pollen hierarchical dominance, expressed
								as relative order, i.e. order/total number of allele*/
	float rel_order_sq;		/*sum of squares of rel_order in order to compute variance*/
	unsigned long age;		/*age of allele with given freq*/
};      

/*defines the structure used to record allelic frequencies distribution	for each class_order
	if pollen hierarchical dominance */
struct type_allele_distrib_class {
	unsigned int count;		/*number of occurence of allele within each class in all simulations*/
	float freqall;			/*  allelic frequency at each class*/
	float freqall_sq;		/*sum of squares of freqall in order to compute variance*/
	unsigned long age;		/*age of allele with given class_order*/
};      

/*defines the structure used to create an allelic tree*/
struct type_tree {   /*as in Hudson 1990*/
	float id;			/*correspond to allele id of terminal nodes and index in tree[index]*/
	long time;
	struct type_tree *desc1;
	struct type_tree *desc2;
	struct type_tree *ancestor;
};

struct typenode{
	long int time;		                /*generation at which node formed*/
	long int alleleid;		                /*id of the allele of node*/
	long int timeallele;		         /*generation at which the allele appeared the first time in the pop.*/
	int mut_up;		                        /*number of mutation events between node and its nearest ancestor in the tree*/
	/*int mut_down; */		                /*number of mutation events between node and all its descendants at the bottom of the tree*/
	int nnext;		                        /*number of nodes below and directly attached to node*/
	/*int nfinal;   */		                        /*number of descendants of node in the last generation at the bottom of the ree*/
	struct typenode *ancestor;	        /*pointer to the nearest ancestor node*/
	struct typenode *descendant[MAX_DESC]; /*array of pointers to each direct descendant nodes*/
	};


/*declarations from file speq.c*/
void build_originalgen(float *newgenzyg,long N2,int model_type,
					   float *allele_order,long *seed);
void build_origgen_monomorph(float *newgenzyg,long N2,int model_type,
					   float *allele_order) ;
int repro(float *oldgenzyg1,float *oldgenzyg2,float *newgenzyg1,float *newgenzyg2,
		  long N2_old,long N2,int Ndemes,double geneflow,double recombrate,int model_of_reproduction,double ovd_s,
			float *allele_order1,long nalleles1,float *effect_geneflow,long *seed);
void applymutation(float *newgenzyg,long N2,double mutrate,long gen,long *lastalleleid,
					int model_of_reproduction,int model_type,float *allele_order,long *nalleles,long *seed);
long get_order(float *allele_order,long na,float allid);
float pollen_dominant(float allele1,float allele2,float *allele_order,long nalleles);
void applyrecomb(float *newgenzyg1,long N2,double recombrate,long *seed);
void update_allele_order(float *allele_order,long *nalleles_order,struct typeallele *allele,long na);
void stat_alleles(float *newgenzyg,long N2,struct typeallele *allele,long *na, float *ne, float *he);
void build_array_allele(float *newgenzyg,long N2,struct typeallele *allele);
long count_alleles(struct typeallele *allele,long N2);
float comp_ssq_freq(struct typeallele *allele,long N2);
void comp_allele_distrib(struct typeallele *allele,long na,int model,
					float *allele_order,long nalleles_order,struct type_allele_distrib *result,
					struct type_allele_distrib_class *result_class,long N2,int **result_all_distrib);
void comp_distrib_fixed_alleles(struct typeallele *allele,long na,int model,
					float *allele_order,int **fixed_all_distrib,long *fixed_all_freq); 
void comp_gen_fixed_alleles(float *newgenzyg,long N2,long **fixed_gen_distrib);
void comp_mate_avail_moth(float *newgenzyg,long N2,int model_of_reproduction,
			float *allele_order,long nalleles,float *mean_avail,float *var_avail) ;
void comp_hoobs(float *newgenzyg,long N2, float *hoobs);

void comp_mate_avail_fath(float *newgenzyg,long N2,int model_of_reproduction,
			float *allele_order,long nalleles,float *mean_avail,float *var_avail);
void get_demes_limit(long N2,int Ndemes,int d,long *limit_inf,long *limit_sup);
void stat_alleles_perdeme(float *newgenzyg,long N2,int Ndemes,struct typeallele **alleledeme,
					long *na_deme,float *naavg_deme,float *ne_deme,float *he_deme);
void compute_GST(struct typeallele *allele,long na,struct typeallele **alleledeme,
				 long *na_deme,long N2,int Ndemes,float *HT,float *HS,float *GST);
void comp_mate_avail_deme_moth(float *newgenzyg,long N2,int Ndemes,int model_of_reproduction,
			float *allele_order,long nalleles,float *mean_avail,float *var_avail);
void checknewgenzyg(float *newgenzyg,long N2,int flag,long gen);
void testran1(long a,long b,long *seed);
void compute_FST_all(struct typeallele *allele,long na,struct typeallele **alleledeme,
				 long *na_deme,long N2,int Ndemes,char *filename, float *FST_WC,float *FST_HR);
void write_sample_distrib_freq(float *newgenzyg,long size,char *filename);

/*declarations from file spgen.c*/
void build_first_all_new(struct typeallele *allele,long nall,struct type_all **all_new,
						long *nall_new);
void applymutation_all(float *newgenzyg,long N2,double mutrate,long gen,long *lastalleleid,
					struct type_all **all,long *nall,int model_of_reproduction,int model_type,
					float *allele_order,long *nalleles,long *seed);
float get_rel_order(float *allele_order,long na,float allid);
struct type_all *search_anc_allele(struct type_all **all,long nall,float oldalleleid);
void update_all_new(struct typeallele *allele,long naobs,struct type_all **all_old,
				long *nall_old,struct type_all **all_new,long *nall_new,long gen,
				int model_type, int nfirstgen,long nall_life,long *count_life,float **life);
int count_lineages(struct type_all **all,long nall);
void comp_Tc_all(struct type_all **all,long nall,long gen,long *Tc, float *Tc_order);
void comp_Td_all(struct type_all **all,long nall,long gen,float *Td,float *Td_order); 
void comp_Td_common(struct typeallele *allele,struct type_all **all,long nall,
						long gen,float *Td,float *Td_order) ;
int comp_maxmut_all(struct type_all **all,long nall);
void comp_Ti_from_all(struct type_all **all,long nall,int na_sample,long gen,int countsimul,
						long **Ti,long **Ti_order,long *seed);
void comp_Ti_from_all_old(struct type_all **all,long nall,int na_sample,long gen,int countsimul,
						long **Ti,long *seed);
void comp_Ti_common(struct typeallele *allele,struct type_all **all,long nall,int na_sample,long gen,int countsimul,
						long **Ti,long *seed);
long get_freqall(struct typeallele *allele,long naobs,float alleleid);
void comp_Tci_cat(struct type_all **all,int nall,long *Tci,int *inode);
void shellsort_int(unsigned long n, int a[]);
void shellsort_long(unsigned long n, long a[]);
void resample_shuffle(long *index, long first,long last, long *seed);
void resample_init(long *index,long first,long last);
void comp_allele_age(struct typeallele *allele,long na,struct type_all **all,long nall,long gen,int model_type,
					float *allele_order,long nalleles_order,
					struct type_allele_distrib *result,struct type_allele_distrib_class *result_class);
void record_order_in_all(float *allele_order,long na,struct type_all **all,long gen);
void comp_Td_hierarchy(struct type_all **all,long nall,long gen,float *allele_order,long **Td_hierarchy); 
void comp_Ti_recurs(struct type_all **all,long nall,long lastT,long *Ti,float *Ti_order);
long get_Tipair(struct type_all *a,struct type_all *b,float *order_par);
void sort_all(struct type_all **all,long nall);
int alli_higher_than_a(struct type_all *alli,struct type_all *al);
void comp_alpha(struct type_all **all,long nall,long gen,float *alpha,long *seed);
void check_all_older_thanTc(struct type_all **all,long nall,float *allele_order,long Tc,long gen,
							int model_of_reproduction,int *flag_order,int *flag_recessive);
void comp_difforder_firstcoal(struct type_all **all,long nall,float *allele_order,float *difforder);
void make_tree(struct type_all **all,long nall,struct type_tree *tree,long gen);
void make_tree_recurs(struct type_all **all,long nall,struct type_tree *tree,long nalltot);
struct type_tree *get_desc_tree(float allid,struct type_tree *tree,long nalltot);
void check_Ti_from_tree(struct type_all **all,long nall,long gen,struct type_tree *tree);
void comp_Ti_from_tree(struct type_tree *tree,long nall,long *Ti);
void comp_Uyenoyama(struct type_tree *tree,long nall,long *T,long *D,float *P,long *S,float *B,
					float *RPT,float *RST,float *RSD,float *RBD);
void comp_Uyenoyama_Td(struct type_tree *tree,long nall,long *T,long *D,float *P,long *S,float *B,
					float *RPT,float *RST,float *RSD,float *RBD);
void store_life_allele(struct type_all *all,long gen,int model_type,
					   long nall_life,long *count_life,float **life);
void compute_dynamics(char* outputfile,long N2,int model_type,long nall,float **life,
					  int nclas,float *freq);
void indexx(long n,float arr[],long indx[]);
void comp_T_deme(int Ndemes,struct typeallele **alleledeme,long *na_deme,long gen,
				struct type_all **all,long nall,struct type_all **alltemp,
				float *record_Tc_deme,float *record_Td_deme,int count_Tc);
void comp_Uyenoyama_deme(int Ndemes,struct typeallele **alleledeme,long *na_deme,long gen,
				struct type_all **all,long nall,struct type_all **alltemp,
				float *RPT_U,float *RST_U,float *RSD_U,float *RBD_U);



/*declarations from file spinout.c*/
void get_initial_conditions(const char *inputfile,int *nseqAl,int *nseqBl,int *nloc, int *nsl, long *ndatasets, char *datafile);
void get_initial_conditions_dynamics(const char *inputfile,int *model_of_reproduction,
								long *N,int *Ndemes, double *geneflow,
								double *mutrate1,double *mutrate2,
								double *recombrate,double *ovd_s,
								int *nsimul,int *na_sample,
								long *nall_life,int *nclas_life,float *freq_life);
void print_initial_conditions(const char *inputfile,int *nseqAl,int *nseqBl,int nloc, int *nsl, long ndatasets, char *datafile);
void write_initial_conditions(char *filename, const char *inputfile,int *nseqAl,int *nseqBl,int nloc, int *nsl, long ndatasets, char *datafile);
void get_dataset(FILE *fp,long dataset,int *nseqAl,int *nseqBl,int nloc,int *nsl,int nsmax,char ***seqAlhs,char ***seqBlhs,char **seqOls, int *nspolyl);
void get_dataset_bug(FILE *fp,long dataset,int *nseqAl,int *nseqBl,int nloc,int *nsl,int nsmax,char ***seqAlhs,char ***seqBlhs,char **seqOls, int *nspolyl);
void write_dataset_fasta(char *filename,int nseqA,int nseqB,int loc, int nsl,char ***seqAlhs,char ***seqBlhs,char **seqOls); 
void write_polyl(char *outputfilename,long dataset,int nloc,struct result_poly *r);
void initialize_write_polyl(char *outputfilename,int nloc);
void write_ABCstat(char *filename,long dataset,int nloc,struct result_poly *r);
void initialize_write_ABCstat(char *filename);

void write_stat_alleles(const char *filename,long gen,long na,float ne,float he);
void writeTi(char *outputfile,int na,int nsimul,long **Ti);
void writeTi_order(char *outputfile,int na,int nsimul,long **Ti_order);
void write_all(char *filename,struct type_all **all,long nall);
void write_life_allele(struct type_all *all,long gen,int model_of_reproduction);
void write_allelic_freq(const char *filename,long gen,long na,struct typeallele *allele,
					float *allele_order, int model);
void avestd_long_miss(long int data[], unsigned long int n, int missing, float *ave, float *std,long *count_obs);
void write_neobs_int(const char *filename,long gen,float ne1,float ne2, long int_equilib);
		
/*declarations from file spexp.c*/
void comp_expectations_gsi(long N, double v,float *na, float *ne,
		float *he,float *acomp, float *limit, float *fsi);
double phix(float x);
double funcj(float x);
double funcj_ln(float x);
void comp_expectations_neutral(long N, double v,float *na, float *ne,float *he);
double phix_neutral(float x);
void iterate_j(long N,double v,float j_in,float *j_out);
double comp_1F1(float a, float b, float s);
double comp_an(float a, long n);
double comp_ln_an(float a, long n);
double funcj_1F1(float j) ;



/*declarations from file spalloc.c*/
struct result_poly *alloc_result_poly(long size);
struct typeallele *alloc_allele(long size);
struct typeallele **alloc_matrix_allele(long nrow,long ncol);
struct type_all **alloc_all(long size);
void alloc_nodes_to_all(struct type_all **all,long size);
void push_node_stack(struct type_all * node);
struct type_all *pop_node_stack(void);
struct type_allele_distrib *alloc_allele_distrib(long size);
struct type_allele_distrib_class *alloc_allele_distrib_class(long size);
struct type_tree *alloc_tree(long size);

int *alloc_int_vector(long size);
long *alloc_longint_vector(long size); 
int *ivector(long nl,long nh);
void free_ivector(int *v, long nl, long nh);
long **longmatrix(long nrl, long nrh, long ncl, long nch);
void free_longmatrix(long **m, long nrl, long nrh, long ncl, long nch); 
float **floatmatrix(long nrl, long nrh, long ncl, long nch);
void free_floatmatrix(float **m, long nrl, long nrh, long ncl, long nch);
int **intmatrix(long nrl, long nrh, long ncl, long nch);
void free_intmatrix( int **m, long nrl, long nrh, long ncl, long nch);
long *lvector(long nl, long nh);
void free_lvector(long  *v, long nl, long nh);
float *vector(long int nl, long int nh);
void free_vector(float *v, long int nl, long int nh);
double *dvector(long int nl, long int nh);
void free_dvector(double *v, long int nl, long int nh);
char ***char3tensor(long nrl, long nrh,long ncl, long nch,long ndl,long ndh);
void free_char3tensor(int ***t,long nrl, long nrh,long ncl, long nch,long ndl,long ndh);
int ***i3tensor(long nrl, long nrh,long ncl, long nch,long ndl,long ndh);
void free_i3tensor(int ***t,long nrl, long nrh,long ncl, long nch,long ndl,long ndh);
char **charmatrix(long nrl, long nrh, long ncl, long nch);
void free_charmatrix( char **m, long nrl, long nrh, long ncl, long nch);


/*declarations from file sputil.c*/
void get_time_short(char *smess);
void get_time_full(char *smess);
void write_info(char *filename,char *info);
void write(char *filename,char *info);
void write_onlytofile(char *filename,char *info);
double sqr(float x);
int roundxa(float x);
long roundxa_long(float x);
void avestdev_longint(long int data[], unsigned long int n,float *ave,float *stdev);
void avestdev_float(float data[], unsigned long int n,float *ave,float *stdev);
double comp_std(long n,double sy,double sy2);
void avestdev_unsignedshort(unsigned short data[], unsigned long int n,float *ave,float *stdev);
void avestd_float_missing(float data[], unsigned long n, float missing, float *ave, float *std);

/*declarations from file sidna.c*/
void compute_polyl(int nseqA,int nseqB,int nsl,int nspolyl,
				   char **seqAlhs,char **seqBlhs,char *seqOls,struct result_poly *r);
double WaldWolfowitz(int number_of_elements, int *IOseq, int std);
void compute_DTajima(int nseq,int nsegsites,float avgpairdif,float *D);

void create_anc_seq(int nnodes,int **site,int nsites,long int *seed); 
void apply_mutation_toseq(int nmut,int *seq,int nsites, long int *seed);
void build_sequences(struct typenode **tree,int nnodes,int nsample,int **site,int nsites,double mutrate2,long int *seed); 
void count_segr_sites(int **site,int nsites,int nseq,int *nsegsites);
void comp_avgpair_nucldiff(int **site,int nsites,int nseq,double *avgpairdif,double *avgpairdif_std);
void Tajima_test(int **site,int nsites,int nseq,float *D,int *nsgsites,float *avgpd,float *avgpd_std); 
float pearson_corr_pi(struct result_poly *resl,int nloc);



 



 













