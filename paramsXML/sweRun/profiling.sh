#!/bin/bash

#SBATCH --job-name=SWEDO16
#SBATCH -n 32
#SBATCH -N 2
#SBATCH --time=36:00:00
#SBATCH --mem=8000
#SBATCH --job-name=sourabh
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=s.kotnala@rug.nl

# dir - present directory
# stren - stochastic force strength 
# dirname - run directory name. For example, currently,
# dirname="run_20x40_{stren(i)}Wz".
# call this file from run directory as  'bash run.sh <no. of cores>' 
# For example : 'bash run.sh 4' will run the code with 4 cores   
dir=${PWD}
# stren=(0.000000001 0.00000001 0.0000001 0.00001 0.0001 0.001)
# stren=(0.01 0.1 1.0 10.0 100.0 1000.0)

node=(2 1 1 1 1)
proc=(32 16 8 4 2)
basis=(70 60 50 40 30 20)
it=(12800 6400 3200 1600 800)

for ((j=0; j<6; j++)); do
  for k in ${basis[@]}; do
    for i in ${it[@]}; do
      dirname="solution_${proc[${j}]}_${k}_${i}"
      # make directory
      if [ -d $dirname ]
      then 
	echo "directory ${dirname} already exists"
      else
	mkdir ${dirname}
      fi
      # params.xml and 'det_sol_20x40.mm' will be copied to each sub-directory 
      cp params.xml ${dirname}/
      cp Amesos_Scalapack.xml ${dirname}/
      cp Amesos_Klu.xml ${dirname}/
      cp AztecOO_GMRES.xml ${dirname}/
      # copy perturbation.xml to sub-directory with changed value of strength
      sed -r '/"Max. Stoch subspace Dimension"/ s/[0-9]+/'${k}'/' StochasticParams.xml > $dirname/StochasticParams.xml
      # copy restart.xml(refers to starting solution) to sub-directory 
      # go to run directoryand make a link to executable in build directory.
      cd ${dirname}
      sed -r '/"Stochastic Iterations"/ s/[0-9]+/'${i}'/' StochasticParams.xml > newStochasticParams.xml
      mv newStochasticParams.xml StochasticParams.xml
      if [ ! -f swedo ]
      then
	ln -s ../../SWEDO swedo
      fi
      # start run 
      echo "running for numproc = ${proc[${j}]}, numbasis = ${k} and NumStochIter = ${i}"
      srun -n ${proc[${j}]} -N ${node[${j}]} ./swedo > output.log
      #srun ./swedo
      echo "done"
      cd ${dir}
    done
  done
done
