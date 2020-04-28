#!/bin/sh
#PBS -N poisson
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=48:r662
#PBS -q mei


module load gcc/5.3.0
module load papi/5.4.1

NODE=662
TOL=10

cd $PBS_O_WORKDIR

if [ ! -d ResPoisson\_$NODE ]; then
    mkdir ResPoisson\_$NODE;
else
    rm -r ResPoisson\_$NODE;
    mkdir ResPoisson\_$NODE;
fi

DIR=ResPoisson\_$NODE

for size in {50,100,500,1000,5000,10000,50000,100000}
    do 
        
        for j in {2,4,8,16,24,32,40,48,64}
        do
            echo "Size:$size NºThreads:$j" >> $DIR/$size\_$j.txt
            ./a.out $size $j >> $DIR/$size\_$j.txt
            echo "" >> $DIR/$size\_$j.txt
        done
    done
