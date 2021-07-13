#!/bin/bash

set OMP_NUM_THREADS=1

items=( 0.34 0.38 0.42 0.46 0.50 0.55 0.60 0.69 0.79 0.90
 1.02 1.18 1.35 1.61 1.93 2.31 2.90 3.63 4.75 6.51 10.0
)

if [ ! -e nldft_sda_slit_openmpi.exe ]; then
	mpic++ -O2 nldft_sda_slit_openmpi.cpp -o nldft_sda_slit_openmpi.exe
fi

if [ ! -d results ]; then
	mkdir results
fi

for w in "${items[@]}"; do
	cp temp_parameters.txt parameters.txt
	echo "Pore width = ${w} [nm]"
	H=`awk -v w=${w} "BEGIN {print w+0.34}"`
	echo "Slit width = ${H} [nm]"
	sed -i "s/XXX/${H}/g" parameters.txt
	mpirun -np $1 ./nldft_sda_slit_openmpi.exe
	mv PP0_vs_Vgamma_data.txt ${w}_data.txt
	cp ${w}_data.txt ./results/${w}_data.txt
	rm parameters.txt
done
