#!/bin/bash

#old version
#items=( 0.34 0.38 0.42 0.46 0.50 0.55 0.60 0.69 0.79 0.90
# 1.02 1.18 1.35 1.61 1.93 2.31 2.90 3.63 4.75 6.51 10.0
#)

#if [ ! -e nldft_sda_slit_HP_double.exe ]; then
	#c++ nldft_sda_slit_HP.cpp -o nldft_sda_slit_HP.exe
	c++ nldft_sda_slit_HP_double.cpp -o nldft_sda_slit_HP_double.exe
#fi

if [ ! -d results ]; then
	mkdir results
fi

#for w in "${items[@]}"; do  #old version
#for ((i=30;i<=1000;i+=1)); do  #new version
for ((i=28;i<=200;i+=1)); do  #new version
	w=`echo $i | awk '{printf "%4.2f", $1/100}'`  #new version
	#echo $w
	cp temp_parameters_H2_HP_double.txt parameters.txt
	echo "Pore width = ${w} [nm]"
	H=`awk -v w=${w} "BEGIN {print w+0.34}"`
	echo "Slit width = ${H} [nm]"
	sed -i "s/XXX/${H}/g" parameters.txt
	./nldft_sda_slit_HP_double.exe
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
