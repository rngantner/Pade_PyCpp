#!/bin/bash

# load modules, submit jobs to run demo_cores.py for
# different numbers of cores. For multiple repetitions,
# just execute this file again- demo_cores.py
# automatically updates the previously computed values.

# import correct modules
module load boost
module load python/2.7.2
#module load icc

# boost_python library path
export LD_LIBRARY_PATH=/cluster/apps/boost/1_47_0_nompi/x86_64/gcc_4.1.2/lib64:$LD_LIBRARY_PATH

# loop over number of cores
for i in 1 2 4
do
  echo ${i}
  export OMP_NUM_THREADS=${i}
  bsub -n ${i} -W 08:00 python demo_cores.py
done

# do this after all jobs are done
#python2 evalcores.py

