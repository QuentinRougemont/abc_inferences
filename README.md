# abc_inferences 

# Purpose
whole pipeline to perform [coalescent](http://www.nature.com/nrg/journal/v3/n5/full/nrg795.html) simulations and  [abc](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002803) inferences

This pipeline was used to reconstruct the demographic history of Atlantic salmon *salmo salar*.

I also provide additional script to reproduce major analysis from the paper about "reconstructing the history of Atlantic salmon *salmo salar* using ABC" [available in early view](https://onlinelibrary.wiley.com/doi/epdf/10.1111/evo.13486)


## major steps:

1. compile the coalescent simulator and stats calculator
2. prepare your input file
3. choose the model you want to run and the appropriate prior
4. run the coalescent sims
5. reshape the data
6. run the abc analysis for : 
 (a) model selection 
 (b) robustness computation 
 (c) parameters estimation and goodness of fit
7. (eventually run the pipeline using only barrier loci)
8. (eventually draw neutral enveloppe to perform outlier detection test - I will add this later )

## 1. Programm compilation
depending on your architecture there is more or less efficient way to run the code.
icc compilation will make code runs faster.
command lines for compilation are provided in:

`compile_msnsam.sh` located in `src/msnsam_code` folder

`compile_mscalc.sh` located in `src/mscalc_code` folder

* Alternatively for The coalescent simulator (msnsam) you can type:

```bash
gcc -O2 -o msnsam  msnsam.c  rand1.c streec.c -lm
gcc -o samplestats tajd.c sample_stats.c -lm
gcc -o mean_std stats.c -lm
```

if you have access to an icc architecture you can speeds up the code :

```bash
icc -O3 -xHost -o ../../bin/msnsam msnsam.c rand1.c streec.c -lm
```

* Then for The summary statistics computer (mscalc):

```bash
gcc *.c  -o mscalc -lm
```

if you have access to an `icc` architecture you can speed up the code :

```bash
icc -O3 -xHost *.c -o ../../bin/mscalc -lm
```

## 2. Prepare input(s) file(s)

The data I used for abc analysis of the *Salmo salar* are in a [companion github repository](https://github.com/QuentinRougemont/abc_inferences_data).   
Download the github repository store the data in `01-salmon_data`
and run:

```
./00-scripts/utility_scripts/prepare_abc.sh 01-salmon_data/salmon.data.ped 01-salmon_data/salmon.data.map
```

all the data will be created in different folders.  
Now these data are ready to perform coalescent simulations. You can clone/paste the folders
`00-scripts/ bin/ results/ and src/`
in all the pairwise folder where you want to run simulations.

The original data are available on dryad:

10.5061/dryad.gm367 and doi:10.5061/dryad.sb601

If you want you can dowload them, merge them and filter them in `plink` using the options:

```
--geno 0.05
--mind 0.05
--noweb
```
this will result in a dataset of ~5000 SNPs and 2035 individuals as provided in the [sister github repository](https://github.com/QuentinRougemont/abc_inferences_data)  

### Runing directly on an input exampl  e
if you are lazy preparing the data i also provide a ready to use example in the folder `00-example_input_abc` for the [sister github repository](https://github.com/QuentinRougemont/abc_inferences_data)   
Then you'll have to choose the models and prior associate to each models.

## 3. Choose models and priors

#### Models

Four classical models are implemented:
1. Strict Isolation (SI)
1. Isolation w. migration (IM)
1. Secondary Contact (SC)
1. Ancient migration (AM)

All models were used for the Salmon data. In addition the following models were implemented for some other project I may publish one day:
1. panmixia
1. bottleneck

The four classic models have different alternatives. See Roux et al. 2013 MBE, Roux et al. 2014 JEB, Roux et al. 2016 Plos Biol for more info
These are:
* NhomoMhomo (=same migration rate M and Ne among loci)
* NhomoMhetero (=same Ne among loci but different M among loci)
* NheteroMhomo (=same M but different Ne among loci )
* NheteroMhetero (= different M and different N among loci)

For the SC models there is also the possibility to test for:
* Periodic Secondary Contacts (PSC) when populations undergone two phases of splits and mixture during the divergence process
* bottleneck post divergence in the daughter populations

### Priors

in `ms` (and `msnsam`) all _Ne_ are scaled by 4\*_Nref_\*µ (\*L) where:  
_Nref_ = size of a reference population  
µ = mutation rate  
L = length of the locus  

so that &#952; = 4\*_N1_\*µ*L/4\*_Nref_\*µ\*L

Similarly &Tau; = Tsplit/4\*_Nref_\*µ\*L

and M = 4\*_Nref_\*m
For more details please look at the `msdoc.pdf` manual before running your inferences.

Once you are okay with `ms` coalescent parameters, you may want to modify the priors:  
All priors on N1, N2, Nancestral, Migration rates, Split times, etc can be edited in the `00-scripts/models/model.*.sh` files.

## 4. run the coalescent sims

Coalescent models are found in `00-scripts/models/model.{1..18}.sh` files.  
Each number corresponds to a different model as follows:  
model.1.sh = SI heterogenous Ne  
model.2.sh = SI homogenous Ne  
model.3.sh = AM heterogenous Ne and heterogenous M  
model.4.sh = AM homogenous   Ne and heterogenous M  
model.5.sh = AM heterogenous Ne and homogenous M  
model.6.sh = AM homogenous Ne and homogenous M   
model.7.sh = SC heterogenous Ne and heterogenous M  
model.8.sh = SC homogenous Ne and heterogenous M  
model.9.sh = SC heterogenous Ne and homogenous M  
model.10.sh = IM heterogenous Ne and heterogenous M  
model.11.sh = IM homogenous Ne and heterogenous M  
model.12.sh = IM heterogenous Ne and homogenous M  
model.13.sh = SC homogenous Ne and homogenous M  
model.14.sh = IM homogenous Ne and homogenous M  
model.15.sh = SCb heterogenous Ne and heterogenous M  
model.16.sh = PSC heterogenous Ne and heterogenous M  
model.17.sh = PSCb heterogenous Ne and heterogenous M  
model.18.sh = PAN heterogenous Ne   
model.00.sh = PAN homogenous Ne  


I've run them using `moab scheduler` on a cluster where I could launch several job in parallel.

This is the purpose of the `00-scripts/models/model_job_array*.sh` scripts that you'll have to customize according to your own cluster.

## 5. reshape the data

run the script `./00-scripts/00.reshape.sh`
This will merge all parallel simulation into one prior file and one summary_statistics file with 1milions row (1 row/simulations) 

## 6. run the abc analysis for :

(a) model selection  
(b) robustness computation  
(c) parameters estimation and goodness of fit

dependencies:
* [R software](https://www.r-project.org/)
* `abc` package : more info [here](https://cran.r-project.org/web/packages/abc/index.html)

##### model selection

run `./00-scripts/01.model.selection.sh`

you'll need some abc background (see references list below)  

##### robustness computation

```bash
cd compute_robustness.abc
run ./01-scripts/01_robustess_job.sh args
```

where arg is the name of the model
either si.1 si.3 am.1 am.2 am.3 am.4 im.1 im.2 im.3 im.4 sc.1 sc.2 sc.3 sc.4

this is very time consumming....

#### parameter estimation

```bash
./00-scripts/03.paramestim.im.all.sh
./00-scripts/03.paramestim.sc.all.sh
./00-scripts/03.paramestim.si.all.sh
./00-scripts/03.paramestim.am.all.sh

./00-scripts/06.goodness_fit.sh #for sc only!!!
```

#### run them all:

```bash
./00-scripts/runall_abc_job.sh #for moab only
```

#### other application :

To do

## Other analyses I did  in the Salmon data:

### Treemix analysis

#### Dependencies

`treemix` available [here](https://bitbucket.org/nygcresearch/treemix/wiki/Home) and relevant reference [here](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002967) and [there](https://www.nature.com/articles/ncomms2140).  
`gsnap` available [here](http://research-pub.gene.com/gmap/)  
`samtools` available [here](https://github.com/samtools/samtools)  
The softwares must be into your path.

####  alignment to the salmon genome

 - download the _Salmo salar_ reference genome [here](ftp://ftp-trace.ncbi.nih.gov/genomes/Salmo_salar/Assembled_chromosomes/seq/)  

 store the data in `02-aditional_analisys/00-genome_alignment/01-input`

 - merge all chromosomes into a single fasta file (I used both the masked and unmasked version, it did not changed much of the results)  
run a simple command like:  
 ```bash
 for i in 00-genome_alignment/01-input/*gz 
 do
    gunzip $i ;  
 done
for i in 00-genome_alignment/01-input/*chrssa*fa 
do 
    cat $i >> 00-genome_alignment/01-input/genome_fa ; 
done
```

**run gsnap**

- build the genome  
`./00-additional_analysis/00-scripts/01.gsnap_build.sh`

- align the sequences  
`./02-additional_analysis/00-scripts/02.gsnap_align.sh`

- then remove hardclip and softclip and order the ped files according the markers positions in the genome  

 (I'll document that soon)

 I provide the reference fasta sequence in the [sister repository](https://github.com/QuentinRougemont/abc_inferences_data)


**run Treemix and f3-tests**

as always, the data (ped files) to reproduce th example are located in the [sister repository](https://github.com/QuentinRougemont/abc_inferences_data) in `02-treemix-data`.  
You can download them, uncompress and store them in `02-additional_analysis/02-treemix/`

```bash
#if you have already used treemix then you'll surely already have your own script to run the software.  
#in this case simply prepare the data as follows  

datafolder=02-additional_analyis/02-treemix
plink --file "$datafolder"/salmon.ordered.V2 --noweb --missing --freq --double-id --allow-extra-chr --chr-set 29 --within cluster.dat --out "$datafolder"/plink  
gzip "$datafolder"/plink.frq.strat  
02-additional_analysis/00-scripts/treemix/plink2treemix.py "$datafolder"/plink.frq.strat.gz "$datafolder"/treemix.frq.gz  

#then you'll know what to do.  

#Alternatively to prepare the input from the ped files and run the f3-test use the following command:  
./02-additional_analysis/00-scripts/treemix/prepare_treemix.sh #will prepare the treemix data and perform f3tests

#Then to run treemix use a script like the one here  
./02-additional_analysis/00-scripts/treemix/run_treemix.sh
#Note that this script is very simplified version, I used a more elaborate script to run in parrallel mode.
```
### Spacemix analysis

You can download the data in the [sister repository](https://github.com/QuentinRougemont/abc_inferences_data) in `03-spacemix-data` and store them in `02-additional_analysis/spacemix/01-data`.

```bash
cd 02-additional_analysis/03-spacemix
#creating input file:
./02-additional_analysis/03-spacemix/00-scripts/00.ped_to_matrice.R path/to/input/input.ped #with input ped the name of the ped_file (with the path/)
./02-additional_analysis/03-spacemix/00-scripts/01.make.allele.freq.cov.2.R path/to/input/matrix #with matrix the matix name and the path

#run spacemix (Here I very simply follow guideline from Bradburg's github)
./02-additional_analysis/03-spacemix/00-scripts/02.run.spacemix.R sample.covariance.matrix sample.size x_y.coordinates
#this takes some times and needs to be run on a cluster

./02-additional_analysis/03-spacemix/00-scripts/03.chek.spacemix.model.R arg1 arg2 arg3

#then plot the results
#This is very straightforward if you follow examples from spacemix manual

```

### Relevant references:

Review about abc:

 * [Beaumont MA, Zhang W, Balding DJ. Approximate Bayesian Computation in Population Genetics. Genetics. 2002;162: 2025–2035.](http://www.genetics.org/content/162/4/2025)
 * [Csillery, K., M.G.B. Blum, O.E. Gaggiotti and O. Francois (2010) Approximate Bayesian Computation (ABC) in practice. Trends in Ecology and Evolution, 25, 410-418.](http://www.sciencedirect.com/science/article/pii/S0169534710000662)

Must Reads about Genome-Wide Heterogeneity:

 * [Roux C, Fraïsse C, Romiguier J, Anciaux Y, Galtier N, Bierne N. Shedding Light on the Grey Zone of Speciation along a Continuum of Genomic Divergence. Moritz C, editor. PLOS Biol. 2016;14: e2000234. doi:10.1371/journal.pbio.2000234](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.2000234)

 * [Roux C, Tsagkogeorga G, Bierne N, Galtier N. Crossing the Species Barrier: Genomic Hotspots of Introgression between Two Highly Divergent Ciona intestinalis Species. Mol Biol Evol. 2013;30: 1574–1587. doi:10.1093/molbev/mst066](https://academic.oup.com/mbe/article/30/7/1574/972546/Crossing-the-Species-Barrier-Genomic-Hotspots-of?cited-by=yes&related-urls=yes&legid=molbiolevol;30/7/1574)

Programms:
 * ms : [Hudson RR. Generating samples under a Wright-Fisher neutral model of genetic variation. Bioinforma Oxf Engl. 2002;18: 337–338.](https://academic.oup.com/bioinformatics/article/18/2/337/225783/Generating-samples-under-a-Wright-Fisher-neutral)

 * msnsam : [Ross-Ibarra J, Tenaillon M, Gaut BS. Historical Divergence and Gene Flow in the Genus Zea. Genetics. 2009;181: 1399–1413. doi:10.1534/genetics.108.097238](http://www.genetics.org/content/181/4/1399)

 * abc R package: [Csilléry K, François O, Blum MGB. abc: an R package for approximate Bayesian computation (ABC). Methods Ecol Evol. 2012;3: 475–479. doi:10.1111/j.2041-210X.2011.00179.x](http://onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2011.00179.x/full)

Salmon paper:
 * [Rougemont, Bernatchez The demographic history of Atlantic Salmon (Salmo salar) across its distribution range reconstructed from Approximate Bayesian Computations, Evolution, 2018, https://doi.org/10.1111/evo.13486](https://onlinelibrary.wiley.com/doi/pdf/10.1111/evo.13486)
 
