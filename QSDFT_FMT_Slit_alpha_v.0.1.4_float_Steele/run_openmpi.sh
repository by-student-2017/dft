#!/bin/bash

export OMP_NUM_THREADS=1

#old version
#items=( 0.34 0.38 0.42 0.46 0.50 0.55 0.60 0.69 0.79 0.90
# 1.02 1.18 1.35 1.61 1.93 2.31 2.90 3.63 4.75 6.51 10.0
#)
#
#items=( 0.385 0.432 0.524 0.640 0.753 0.889
# 1.05 1.25 1.48 1.76 2.49 2.97 3.54 4.22 5.04 60.2 72.0 86.0 103.0
#)

#if [ ! -e qsdft_fmt_slit_steele_openmpi.exe ]; then
	mpic++ -O2 qsdft_fmt_slit_steele_openmpi.cpp -o qsdft_fmt_slit_steele_openmpi.exe
#fi

if [ ! -d results ]; then
	mkdir results
fi

run_command="mpirun -np $1 ./qsdft_fmt_slit_steele_openmpi.exe"
#for w in "${items[@]}"; do  #old version
for ((i=42;i<=1000;i+=2)); do  #new version
	w=`echo $i | awk '{printf "%4.2f", $1/100}'`  #new version
	#echo $w
	cp temp_parameters.txt parameters.txt
	echo "Pore width = ${w} [nm]"
	H=`awk -v w=${w} "BEGIN {print w+0.3665}"` #0.34
	echo "Slit width = ${H} [nm]"
	sed -i "s/XXX/${H}/g" parameters.txt
	${run_command}
	if [ -e PP0_vs_Vgamma_data_vs.txt ]; then
	  mv PP0_vs_Vgamma_data_vs.txt ${w}_data_vs.txt
	  mv PP0_vs_Vgamma_data_ls.txt ${w}_data_ls.txt
	fi
	if [ -e Pa_vs_Vgamma_data_vs.txt ]; then
	  mv Pa_vs_Vgamma_data_vs.txt ${w}_data_vs.txt
	  mv Pa_vs_Vgamma_data_ls.txt ${w}_data_ls.txt
	fi
	if [ -e atm_vs_Vgamma_data_vs.txt ]; then
	  mv atm_vs_Vgamma_data_vs.txt ${w}_data_vs.txt
	  mv atm_vs_Vgamma_data_ls.txt ${w}_data_ls.txt
	fi
	cp ${w}_data_vs.txt ./results/${w}_data_vs.txt
	cp ${w}_data_ls.txt ./results/${w}_data_ls.txt
	rm parameters.txt
done
