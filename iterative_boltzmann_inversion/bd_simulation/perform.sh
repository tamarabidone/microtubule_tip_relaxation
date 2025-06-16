#!/bin/bash

python3 prepare_microtubule_gdp_tab.py > microtubule_gdp_tab.txt
python3 prepare_microtubule_gtp_tab.py > microtubule_gtp_tab.txt
lmp_serial -var seed1 `bash -c 'echo $RANDOM'` -var seed2 `bash -c 'echo $RANDOM'` -var seed3 `bash -c 'echo $RANDOM'` -in in.microtubule_gdp_tab.txt
lmp_serial -var seed1 `bash -c 'echo $RANDOM'` -var seed2 `bash -c 'echo $RANDOM'` -var seed3 `bash -c 'echo $RANDOM'` -in in.microtubule_gtp_tab.txt