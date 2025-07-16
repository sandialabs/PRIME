#!/bin/bash

PATH2SRC=$1
regions=(B SF V)

for i in ${regions[@]}; do

  echo "Postprocessing region ${i}"
  froot=`python3 -c "import sys, json; print(\"\".join(json.load(open('setup_"${i}".json'))['regioninfo']['region_tag']))"`
  fmcmc=${froot}"_mcmc.h5"
  if ! test -f ${fmcmc}; then
    echo "${fmcmc} does not exist."
    python ${PATH2SRC}/prime_run.py setup_${i}.json
  fi
  python ${PATH2SRC}/prime_compute_epi_inf_curves.py setup_${i}.json 0

  mv _infection.csv   ${i}_infection.csv
  mv _infection.pdf   ${i}_infection.pdf
  mv _infection.h5    ${i}_infection.h5
  mv _forecast_pp.csv ${i}_forecast.csv
  mv _forecast.pdf    ${i}_forecast.pdf
  mv _forecast.h5     ${i}_forecast.h5

done



