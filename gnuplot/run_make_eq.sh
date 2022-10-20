#!/bin/bash

c++ -O2 make_eq.cpp -o make_eq.exe

search_file="_data_vs.txt"
i=1
for file_name in *${search_file}; do
	echo $file_name" "$i
	w=${file_name%%_*}
	cp ${w}_data_vs.txt data_vs.txt
	cp ${w}_data_ls.txt data_ls.txt
	./make_eq.exe
	mv data_eq.txt ${w}_data_eq.txt
	rm data_vs.txt data_ls.txt
	i=$(( $i+1 ))
done
