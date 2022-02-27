#!/bin/bash

#SBATCH --job-name=SWEDO16
#SBATCH -n 16
#SBATCH -N 1
#SBATCH --time=8:00:00
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
# node and proc are 2-tuples hence same index ('j') for the both.
proc=(4 2)
node=(2 1 1 1 1)
basis=(4 8)
it=(1600 3200)
#slvrPkg and mnslvrPkg are 2-tuple called by indx 'l'.
slpkg=("Amesos2" "Amesos" "Belos")
mnslpkg=("Direct" "Direct" "Iterative")
len=${#slpkg[@]}
plen=${#proc[@]}

for ((j=0; j<${plen}; j++)); do
  for k in ${basis[@]}; do
    for i in ${it[@]}; do
      for ((l=0; l<${len}; l++)); do
	dirname="solution_${proc[${j}]}_${k}_${i}_${mnslpkg[$l]}_${slpkg[$l]}"
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
	sed -r '/"Solver Package"/ s/Belos/'${slpkg[$l]}'/' newStochasticParams.xml > StochasticParams.xml
	#mv StochasticParams.xml StochasticParams.xml
	sed -r '/"Solver Type"/ s/Direct/'${mnslpkg[$l]}'/' params.xml > newparams.xml
	mv newparams.xml params.xml
	if [ ! -f swedo ]
	then
	  ln -s ../../SWEDO swedo
	fi
	# start run 
	echo "running for numproc = ${proc[${j}]}, numbasis = ${k} and NumStochIter = ${i}"
	srun -n ${proc[${j}]} -N ${node[${j}]} ./swedo > output.log
	#mpirun -np ${proc[${j}]} ./swedo
	#srun ./swedo
	echo "done"
	cd ${dir}
      done
    done
  done
done
