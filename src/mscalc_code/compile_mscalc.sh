#module load compilers/gcc/4.8
#gcc *.c  -o mscalc -lm
#mv mscalc ../../bin

icc -O3 -xHost *.c -o ../../bin/mscalc -lm
