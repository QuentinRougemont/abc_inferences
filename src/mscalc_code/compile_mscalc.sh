#module load compilers/gcc/4.8
#gcc *.c  -o mscalc -lm
#mv mscalc ../../bin

#module load  compilers/intel/14.0
icc -O3 -xHost *.c -o ../../bin/mscalc -lm
