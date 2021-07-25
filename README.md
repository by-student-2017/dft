# dft (CC0 license)


## Development version.
- For the first time in the world, a significant level of NLDFT code was built in C++ with a CC0 license. 
- I want someone to develop QSDFT, equilibrium, NLDFT(FMT), 2D-NLDFT and RSDFT code with CC0 license and publish them.


##  NLDFT (SDA) beta version 0.9.3 (July/13/2021)
- add "grand potential" calculation routine and adsorption amount in other units (alpha version)
- The equilibrium point is obtained from the intersection of the ground potentials during adsorption and desorption.


##  NLDFT (FMT) beta version 0.9.3 (July/20/2021)
	sudo apt update
	sudo apt -y g++
	cd NLDFT_FMT_Slit_v0.9.3
	gedit temp_parameters.txt
	chmod +x run.sh
	./run.sh
- With only two weeks to work on it, there's still a lot of room to optimize and improve the code, but we hope that building on this achievement will improve the research environment for many people.


##  NLDFT (SDA) beta version 0.9.2 (OpenMPI parallel calculation) (July/10/2021)
	sudo apt update
	sudo apt -y openmpi-bin
	cd NLDFT_SDA_Slit_v0.9.2
	gedit temp_parameters.txt
	chmod +x run_openmpi.sh
	./run_openmpi.sh 4
- e.g., 4 CPU calculation. please, # of CPU is <= (nstep-1)/2 and <= ndmesh = nrmesh.


##  NLDFT (SDA) beta version 0.9.2 (OpenMP parallel calculation) (July/9/2021)
	sudo apt update
	sudo apt -y g++
	cd NLDFT_SDA_Slit_v0.9.2
	gedit temp_parameters.txt
	chmod +x run_openmp.sh
	./run_openmp.sh 2
- e.g., 2 CPU calculation. please, # of CPU is <= (nstep-1)/2 and <= ndmesh = nrmesh. Most "for" loops could not be parallelized.


##  NLDFT (SDA) beta version 0.9.0 (July/7/2021)
	sudo apt update
	sudo apt -y g++
	cd NLDFT_SDA_Slit_v0.9.0
	gedit temp_parameters.txt
	chmod +x run.sh
	./run.sh
- With only six days to work on it, there's still a lot of room to optimize and improve the code, but we hope that building on this achievement will improve the research environment for many people.


## NLDFT (SDA) alpha version 0.8.0 (July/6/2021)
	sudo apt update
	sudo apt -y g++
	cd NLDFT_SDA_Slit_alpha_v0.8.0
	c++ -O2 nldft.cpp -o nldft.exe
	./nldft.exe
	nldft routine is not good !
	Maxwell constraction routine is not good !


## QSDFT (FMT) alpha version 0.1.0 (July/12/2021)
	sudo apt update
	sudo apt -y g++
	cd QSDFT_FMT_Slit_alpha_v.0.1.0
	c++ -O2 qsdft.cpp -o qsdft.exe
	./qsdft.exe
	qsdft routine is not good !
	Maxwell constraction routine is not good !
- This QSDFT just worked. The operation has not been examined in detail yet.
- I'll leave it here as it may be an important step.
- The parameters of non-graphited carbon black BP-280 were used.