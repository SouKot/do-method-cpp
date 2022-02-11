#!/bin/bash

#SBATCH --job-name=SWEDO16
#SBATCH -n 16
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --mem=4000
#SBATCH --job-name=
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=

# dir - present directory
# stren - stochastic force strength 
# dirname - run directory name. For example, currently,
# dirname="run_20x40_{stren(i)}Wz".
# call this file from run directory as  'bash run.sh <no. of cores>' 
# For example : 'bash run.sh 4' will run the code with 4 cores   
dir=${PWD}
# stren=(0.000000001 0.00000001 0.0000001 0.00001 0.0001 0.001)
# stren=(0.01 0.1 1.0 10.0 100.0 1000.0)
#proc=(4 2)
#basis=(70 60 50 40 30 20 10)
#it=(12800 6400 3200 1600 800)
node=(1 1 1 1 1) # number of nodes 
proc=(16 8 4 2 1) # number of total processors
basis=(40) # number of stochastic basis
it=(12800) # number of stochastic iterations
nx=(12800) # number of grid points

for j in ${proc[@]}; do
  for k in ${basis[@]}; do
    for i in ${it[@]}; do
      dirname="solution_proc_${j}_m_${k}_iter_${i}_nx_${nx[0]}"
      # make directory
      if [ -d $dirname ]
      then 
	echo "directory ${dirname} already exists"
      else
	mkdir ${dirname}
      fi
      # params.xml and 'det_sol_20x40.mm' will be copied to each sub-directory 
      cp *.xml ${dirname}/
      # copy perturbation.xml to sub-directory with changed value of strength
      sed -r '/"Max. Stoch subspace Dimension"/ s/[0-9]+/'${k}'/' StochasticParams.xml > $dirname/StochasticParams.xml
      # copy restart.xml(refers to starting solution) to sub-directory 
      # go to run directoryand make a link to executable in build directory.
      cd ${dirname}
      sed -r '/"Stochastic Iterations"/ s/[0-9]+/'${i}'/' StochasticParams.xml > newStochasticParams.xml
      mv newStochasticParams.xml StochasticParams.xml
      sed -r '/"nx"/ s/[0-9]+/'${nx[0]}'/' params.xml > newparams.xml
      mv newparams.xml params.xml
      if [ ! -f brgrdo ]
      then
	ln -s ../../BRGRDO brgrdo
      fi
      # start run 
      echo "running for numproc = ${j}, numbasis = ${k}, nx = ${nx[0]} and NumStochIter = ${i}"
      srun -n ${proc[${j}]} -N ${node[${j}]} ./brgrdo > output.log
      #mpirun -np ${j} brgrdo 2>&1 | tee output.log
      echo "done"
      echo "*********************************************************************************"
      cd ${dir}
    done
  done
done
