#module load compilers/gcc/4.8
#gcc -O2 -o msnsam  msnsam.c  rand1.c streec.c -lm
#gcc -o samplestats tajd.c sample_stats.c -lm
#gcc -o mean_std stats.c -lm
#mv msnsam ../../bin

icc -O3 -xHost -o ../../bin/msnsam msnsam.c rand1.c streec.c -lm
