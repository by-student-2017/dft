# dft (CC0 license)

## Development version.
- For the first time in the world, a significant level of NLDFT code was built in C++ with a CC0 license. 
- I want someone to develop QSDFT and spinodal versions by CC0 license.

## beta (OpenMP parallel calculation) version (July/9/2021)
	cd NLDFT_SDA_Slit_v1.1.1
	gedit temp_parameters.txt
	chmod +x run_openmp.sh
	./run_openmp.sh 2
- e.g., 2 CPU calculation. please, # of CPU is <= (nstep-1)/2. Most "for" loops could not be parallelized.

## beta version (July/7/2021)
	cd NLDFT_SDA_Slit_v1.0.0
	gedit temp_parameters.txt
	chmod +x run.sh
	./run.sh
- With only six days to work on it, there's still a lot of room to optimize and improve the code, but we hope that building on this achievement will improve the research environment for many people.


## alpha version (July/6/2021)
	c++ -O2 nldft.cpp -o nldft.exe
	./nldft.exe
	nldft routine is not good !
	Maxwell constraction routine is not good !