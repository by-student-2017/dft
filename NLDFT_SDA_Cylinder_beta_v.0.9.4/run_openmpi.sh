#!/bin/bash

OMP_NUM_THREADS=1

#items=( 0.34 0.38 0.42 0.46 0.50 0.55 0.60 0.69 0.79 0.90
# 1.02 1.18 1.35 1.61 1.93 2.31 2.90 3.63 4.75 6.51 10.0
#)

if [ ! -e nldft_sda_cylinder_openmpi.exe ]; then
	mpic++ -O2 nldft_sda_cylinder_openmpi.cpp -o nldft_sda_cylinder_openmpi.exe
fi

if [ ! -d results ]; then
	mkdir results
fi

run_command="mpirun -np $1 ./nldft_sda_cylinder_openmpi.exe"
#for w in "${items[@]}"; do
for ((i=34;i<=1000;i+=2))
do
	w=`echo $i | awk '{printf "%4.2f", $1/100}'`
	#echo $w
	cp temp_parameters.txt parameters.txt
	echo "Pore width = ${w} [nm]"
	H=`awk -v w=${w} "BEGIN {print w+0.34}"`
	echo "Slit width = ${H} [nm]"
	sed -i "s/XXX/${H}/g" parameters.txt
	${run_command}
	mv PP0_vs_Vgamma_data_vs.txt ${w}_data_vs.txt
	mv PP0_vs_Vgamma_data_ls.txt ${w}_data_ls.txt
	cp ${w}_data_vs.txt ./results/${w}_data_vs.txt
	cp ${w}_data_ls.txt ./results/${w}_data_ls.txt
	rm parameters.txt
done
