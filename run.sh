#!/bin/sh
#PBS -N AP
#PBS -l walltime=0:15:00
#PBS -l nodes=1:ppn=48:r662


RUNS=100

module load gcc/5.3.0
module load gcc/7.2.0

cd $PBS_O_WORKDIR

make

for CITIES in {512,1024,2048,4096,8192,16384}
do
  if [ $CITIES = 16384]; then
    ./TSP $CITIES 10
  else
    ./TSP $CITIES $RUNS
  fi
done
