#!/bin/bash

# old version
#items=( 0.34 0.38 0.42 0.46 0.50 0.55 0.60 0.69 0.79 0.90
# 1.02 1.18 1.35 1.61 1.93 2.31 2.90 3.63 4.75 6.51 10.0 
#)

search_file="_data_vs.txt"

set ra
set ra2
set nm

touch temp.txt
#echo -n > temp.txt # another method

set i
ra="1.1,1.2"
nm=""
#for w in "${items[@]}"; do # old version
for file_name in *${search_file}; do # new version
	w=${file_name%%_*} # new version
	sed -e '1,2d' ${file_name} > temp.txt
	i=$(( $i+1 ))
	if [ $i == 1 ]; then
	  cp temp.txt kernel.csv
	elif [ $i == 2 ]; then
	  join -o 1.1,1.3,2.3 kernel.csv temp.txt > kernel_temp.csv
	else
	  ra=${ra}",1."${i}
	  ra2=${ra}",2.3"
	  join -o ${ra2} kernel.csv temp.txt > kernel_temp.csv
	  mv kernel_temp.csv kernel.csv
	fi
	nm=${nm}","${w}
done
sed -i "1s/^/${nm}\n/" kernel.csv
rm temp.txt
