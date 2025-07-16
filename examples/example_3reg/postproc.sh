#!/bin/bash

PATH2SRC=$1
regions=(B SF V)
regcat=`(IFS=; echo "${regions[*]}")`


echo "Postprocessing regions: ${regions[@]}"
froot=`python3 -c "import sys, json; print(\"\".join(json.load(open('setup_"${regcat}".json'))['regioninfo']['region_tag']))"`
fmcmc=${froot}"_mcmc.h5"
if ! test -f ${fmcmc}; then
  echo "${fmcmc} does not exist."
  python ${PATH2SRC}/prime_run.py setup_${regcat}.json
fi

regnames=(`python3 -c "import sys, json; print(\" \".join(json.load(open('setup_"${regcat}".json'))['regioninfo']['region_tag']))"`)
for i in "${!regions[@]}"; do 
  echo "Postprocessing ${regnames[$i]}"
  ln -fs ${regnames[$i]}_plus.dat plus10.dat
  python ${PATH2SRC}/prime_compute_epi_inf_curves.py setup_${regcat}.json ${i}
  mv _infection.csv   ${regions[$i]}_infection.csv
  mv _infection.pdf   ${regions[$i]}_infection.pdf
  mv _infection.h5    ${regions[$i]}_infection.h5
  mv _forecast_pp.csv ${regions[$i]}_forecast.csv
  mv _forecast.pdf    ${regions[$i]}_forecast.pdf
  mv _forecast.h5     ${regions[$i]}_forecast.h5
done
