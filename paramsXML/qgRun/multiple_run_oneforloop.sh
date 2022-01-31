#!/bin/bash

# dir - present directory
# stren - stochastic force strength 
# dirname - run directory name. For example, currently,
# dirname="run_20x40_{stren(i)}Wz".
# call this file from run directory as  'bash run.sh <no. of cores>' 
# For example : 'bash run.sh 4' will run the code with 4 cores   
dir=${PWD}

#re=(35.0 32.0 33.0 34.0)
#re=(32.0 30.45 29.86)
#re=(33.0 32.8 32.5 32.0 31.5)
re=(40.0)
m=(6 4 4 5 6 6 6)
sf=(1.0 1.0 1.0 1.0 1.0 0.0)
topo=(0.0 0.0 0.0 0.0 0.0 0.0)
initT=( 0.0 10.0 0.0 0.0)
finalT=(10.0 20.0 20.0)
dt=(0.001 0.001 0.005 20.0)
isstable=(1 0 0 0 0 1 0 0 0)
#initMeanSolFile=("None" "None" "None" "None" "None" "None")
initMeanSolFile=("None" "MeanSolnFinal.mm" "MeanSolnFinal.mm" "MeanSolnFinal.mm" )
#initBasisFile=("None" "None" "None" "None"  "None" "None" "None" "None" )
initBasisFile=("None" "V_BASE.mm" "V_BASE.mm" "V_BASE.mm")
## NOTE: You should first change "Type of Coeff. file." to what you want (i.e. if it
## is "None", "CoeffMatrix" or "Variance" file that your are supplying) in
## StochasticParams.xml file. Then you supply the file name below accordingly!!!
#coefFlType=( "None" "None" "None" "Variance" "None" "None" "None" "None" "Variance")
#coefFlType=("Variance" "Variance" "Variance" "Variance" "Variance" "Variance" "Variance")
coefFlType=("None" "CoeffMatrix" "CoeffMatrix" "CoeffMatrix" "CoeffMatrix" "CoeffMatrix")
#initCoefFile=("initvar.mm" "initvar.mm" "initvar.mm" "initvar.mm" "initvar.mm" "initvar.mm" )
initCoefFile=("yT_3.17999992e+01.mm" "YTrans_COEFF.mm" "YTrans_COEFF.mm" "YTrans_COEFF.mm")
len=${#re[@]}
for ((j=0; j<len; j++)); do
  echo "running for re = ${re[$j]}, numbasis = ${m[$j]}, StchFrcStren = ${sf[${j}]} and topo_coef = ${topo[$j]}"
  dirname="re_${re[$j]}_m_${m[$j]}_sf_${sf[${j}]}_topo_${topo[$j]}_isstable_${isstable[$j]}_multistable"
  # make directory
  if [ -d $dirname ]
  then 
    echo "directory ${dirname} already exists"
  else
    echo ${dirname}
    mkdir ${dirname}
  fi
  # params.xml and 'det_sol_20x40.mm' will be copied to each sub-directory 
  cp *.xml ${dirname}/
  # copy perturbation.xml to sub-directory with changed value of strength
  sed -r '/"Max. Stoch subspace Dimension"/ s/[0-9]+/'${m[$j]}'/' StochasticParams.xml > $dirname/StochasticParams.xml
  # copy restart.xml(refers to starting solution) to sub-directory 
  # go to run directoryand make a link to executable in build directory.
  cd ${dirname}
  sed -r '/"StochFrc Strength"/ s/[0-9]\.[0-9]+/'${sf[${j}]}'/' StochasticParams.xml > newStochasticParams.xml
  sed -r '/"StochBasisFile"/ s/None/'${initBasisFile[$j]}'/' newStochasticParams.xml > NewStochasticParams.xml
  sed -r '/"StochCoefFile"/ s/None/'${initCoefFile[$j]}'/' NewStochasticParams.xml > newStochasticParams.xml
  sed -r '/"Type of Coeff File"/ s/None/'${coefFlType[$j]}'/' newStochasticParams.xml > NewStochasticParams.xml
  mv NewStochasticParams.xml StochasticParams.xml
  sed -r '/"Reynolds Number"/ s/[0-9]\.[0-9]+/'${re[$j]}'/' params.xml > newparams.xml
  sed -r '/"Topography"/ s/[0-9]\.[0-9]+/'${topo[$j]}'/' newparams.xml > Newparams.xml
  sed -r '/"End Time"/ s/[0-9]\.[0-9]+/'${finalT[$j]}'/' Newparams.xml > newparams.xml
  sed -r '/"Start Time"/ s/[0-9]\.[0-9]+/'${initT[$j]}'/' newparams.xml > Newparams.xml
  sed -r '/"Stable Solution"/ s/[0-9]/'${isstable[$j]}'/' Newparams.xml > newparams.xml
  sed -r '/"Initial Solution File"/ s/None/'${initMeanSolFile[$j]}'/' newparams.xml > Newparams.xml
  sed -r '/"Step Size"/ s/[0-9]\.[0-9]+/'${dt[$j]}'/' Newparams.xml > newparams.xml
  mv newparams.xml params.xml
  if [ ! -f qgdo ]
  then
    ln -s ../../QGDO qgdo
  fi
  # start run 
  echo "running for re = ${re[$j]}, numbasis = ${m[$j]}, StchFrcStren = ${sf[${j}]} and topo_coef = ${topo[$j]}"
  strt=$SECONDS
  ./qgdo
  end=$SECONDS
  duration=$(( end - strt ))
  echo "It took $duration seconds to complete run"

  matcode="analyzeNplot('MeanSolnFinal.mm','V_BASE.mm','YTrans_COEFF.mm',64,64,20,10000,1:${m[$j]},'T_${finalT}');\
         cd ..;read_timeseries_Eyy_N_plot('${dirname}',${m[$j]});exit"
  #matcode="cd ..;read_timeseries_Eyy_N_plot('${dirname}',${m[$j]});exit"
  matlab -nodesktop -r "${matcode}"
  #srun ./swedo
  echo "done"
  cd ${dir}
done
