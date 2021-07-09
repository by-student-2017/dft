# dft (CC0 license)


## Development version.
- For the first time in the world, a significant level of NLDFT code was built in C++ with a CC0 license. 
- I want someone to develop QSDFT and spinodal versions by CC0 license.


## beta (OpenMPI parallel calculation) version 0.9.2 (July/10/2021)
	sudo apt update
	sudo apt -y openmpi-bin
	cd NLDFT_SDA_Slit_v1.1.1
	gedit temp_parameters.txt
	chmod +x run_openmpi.sh
	./run_openmpi.sh 4
- e.g., 4 CPU calculation. please, # of CPU is <= (nstep-1)/2 and <= ndmesh = nrmesh.


## beta (OpenMP parallel calculation) version 0.9.2 (July/9/2021)
	sudo apt update
	sudo apt -y g++
	cd NLDFT_SDA_Slit_v1.1.1
	gedit temp_parameters.txt
	chmod +x run_openmp.sh
	./run_openmp.sh 2
- e.g., 2 CPU calculation. please, # of CPU is <= (nstep-1)/2 and <= ndmesh = nrmesh. Most "for" loops could not be parallelized.


## beta version 0.9.0 (July/7/2021)
	sudo apt update
	sudo apt -y g++
	cd NLDFT_SDA_Slit_v1.0.0
	gedit temp_parameters.txt
	chmod +x run.sh
	./run.sh
- With only six days to work on it, there's still a lot of room to optimize and improve the code, but we hope that building on this achievement will improve the research environment for many people.


## alpha version (July/6/2021)
	sudo apt update
	sudo apt -y g++
	c++ -O2 nldft.cpp -o nldft.exe
	./nldft.exe
	nldft routine is not good !
	Maxwell constraction routine is not good !