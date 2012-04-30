#!/bin/bash

# execute demo_cores for different numbers of cores
# and the following number of repetitions:
reps=5

# loop over cores
for i in 1 2 4
do
  echo ${i}
  for j in `seq $reps`
  do
    #export OMP_NUM_THREADS=${i}
    OMP_NUM_THREADS=${i} python2 demo_cores.py
  done
done

# do this after all jobs are done
python2 evalcores.py

