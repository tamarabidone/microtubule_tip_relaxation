#!/bin/bash

python3 dist_analysis.py gdp
python3 dist_analysis.py gtp
python3 angle_analysis.py gdp
python3 angle_analysis.py gtp
python3 dihe_analysis.py gdp
python3 dihe_analysis.py gtp
python3 compute_histograms.py
python3 plot_histograms.py
python3 compute_pot.py $2
python3 plot_pot.py
python3 compute_cylindrical.py gdp
python3 compute_cylindrical.py gtp
python3 plot_cylindrical.py
python3 prepare_potentials.py $1