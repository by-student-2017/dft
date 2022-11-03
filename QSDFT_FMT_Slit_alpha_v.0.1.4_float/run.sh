#!/bin/bash

items=( 0.385 0.432 0.524 0.640 0.753 0.889
 1.05 1.25 1.48 1.76 2.49 2.97 3.54 4.22 5.04 60.2 72.0 86.0 103.0
)

#if [ ! -e qsdft_fmt_slit.exe ]; then
	c++ -O2 qsdft_fmt_slit.cpp -o qsdft_fmt_slit.exe
#fi

if [ ! -d results ]; then
	mkdir results
fi

for w in "${items[@]}"; do
	cp temp_parameters.txt parameters.txt
	echo "Pore width = ${w} [nm]"
	H=`awk -v w=${w} "BEGIN {print w+1.555}"` #1.555
	echo "Slit width = ${H} [nm]"
	sed -i "s/XXX/${H}/g" parameters.txt
	./qsdft_fmt_slit.exe
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
