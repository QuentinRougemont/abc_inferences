#define MSNSAM  0	/*0 if read msnsam results from file  1 if launches Msnsam*/
#define NARMAX  50	/*max number of arguments to msnsam*/

#define GENEALOGY 0	/*0 if no allelic genealogy to perform, 1 if genealogy to perform*/
#define POLLEN_PACKET_MIGRATION 0	/*1 if pollen packet migration model, 0 if single pollen migration model*/
#define START_TRACK	50000			/*20000 number of generations to reach equilibrium*/
#define END_TRACK	80000			/*50000number of generations to reach equilibrium*/
#define N_STEP_MULTIPLE 50	/*number of replicates per run for recording averages*/
#define INITIAL_LOCUS2_MONOMORPHIC 1	/*1 if locus 2 is monomorphic at the initial generation*/
#define COMMON_FREQ	5	/*minimum number of gene copies for an allele to be considered as common*/
#define MAXTRIALS_MOTHER /*500*/ 500 /*max number of trials of new mother gene in procedure repro*/
#define MAXTRIALS_FATHER /*1*/ 500 	/*max number of trials of new pollen in procedure repro
										1: fecundity selection
										500: normal model*/
#define FIXED_ALLELES	0		/*if 1 fixes the number of alleles to get equilibrium freq*/
#define N_FIXED_ALLELES 3		/*number of alleles in the fixed alleles model*/
#define SAMPLE_COMMON_ALLELES  0 /*if 1 then will only sample alleles with a frequency higher than 
									COMMON_FREQ in order to compute Td and interval coalescence times*/

/*defines which additional type of information will be written*/
#define DEBUG_EQUILIBRIUM  0		/*if 1 writes naobs,neobs and heobs for each generation under equilibrium phase*/
#define DEBUG_INT_EQUILIBRIUM  0	/*if 1 writes the average of neobs over each interval*/
#define RECORD_ALLELE_DISTRIB  1	/*if 1 writes the distribution of allele frequencies over all replicate simulations*/
#define N_CLASS_ORDER 5			/*number of classes of alleles with dominance models when record allelic distribution*/
#define RECORD_GEN	0				/*if 1 writes intermediary genealogical information*/
#define N_STEP_GEN  10			/*number of times will write intermediary genealogical info */
#define WARN_MAXTRIALS_FATHER 0	/*if 1 warns when number of pollen trials is higher than MAXTRIALS_FATHER*/
#define WRITE_LIFE_ALLELE 0		/*if 1 records statistics about frequency and order of individual alleles across generations*/
#define WRITE_INTERM_TOPOLOGY 1	/*if 1 writes allelic topology for each replicate*/
#define RECORD_MATE_AVAILABILITY 0	/*if 1 records mean and variance of mate availability*/
#define RECORD_ALPHA 1	/*if 1 then records values of alpha the substitution rate*/
#define UYENOYAMA 1	/*if 1 then computes Marcy's statistics*/
#define RECORD_MULTIPLE 1	/*if 1 then records averages over multiple replicates per run*/

