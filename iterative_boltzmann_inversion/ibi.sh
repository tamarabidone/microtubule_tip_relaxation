#!/bin/bash

for i in {1..60};
	do echo ROUND "$i"
	cd bd_simulation/
	bash perform.sh
	mkdir iteration_"$i"/
	cp *.xyz iteration_"$i"/
	cd ..
	cd bd_trajectory_analysis
	bash perform.sh "$i" 0.1
	mkdir data/iteration_"$i"/
	mkdir figures/iteration_"$i"/
	# cp -r data/dists/ data/iteration_"$i"/
	# cp -r data/angles/ data/iteration_"$i"/
	# cp -r data/dihes/ data/iteration_"$i"/
	cp -r data/hist/ data/iteration_"$i"/
	cp -r data/pot/ data/iteration_"$i"/
	cp -r data/cylindrical/ data/iteration_"$i"/
	cp -r figures/hist/ figures/iteration_"$i"/
	cp -r figures/pot/ figures/iteration_"$i"/
	cp -r figures/cylindrical/ figures/iteration_"$i"/
	cd ..
done
