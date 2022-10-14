#!/bin/bash

#items=( 0.34 0.42 0.56 0.60 0.75 0.80 0.89 0.99 1.00
# 1.18 1.35 1.61 1.93 2.31 2.90 3.63 4.75 6.51 10.0
#)

if [ ! -e qsdft_fmt_slit.exe ]; then
	//c++ -O2 qsdft_fmt_slit.cpp -o qsdft_fmt_slit.exe
	c++ -O2 qsdft_fmt_slit_steele -o qsdft_fmt_slit.exe
fi

if [ ! -d results ]; then
	mkdir results
fi

#for w in "${items[@]}"; do
for ((i=36;i<=1000;i+=2)); do
	w=`echo $i | awk '{printf "%4.2f", $1/100}'`
	#echo $w
	cp temp_parameters.txt parameters.txt
	echo "Pore width = ${w} [nm]"
	H=`awk -v w=${w} "BEGIN {print w+0.195}"`
	echo "Slit width = ${H} [nm]"
	sed -i "s/XXX/${H}/g" parameters.txt
	./qsdft_fmt_slit.exe
	mv PP0_vs_Vgamma_data_vs.txt ${w}_data_vs.txt
	mv PP0_vs_Vgamma_data_ls.txt ${w}_data_ls.txt
	cp ${w}_data_vs.txt ./results/${w}_data_vs.txt
	cp ${w}_data_ls.txt ./results/${w}_data_ls.txt
	rm parameters.txt
done
