# make
#nldft: nldft.o maxwell_construction.o
#	c++ -Wall -o nldft.o maxwell_construction.o
nldft: nldft.o
	c++ -Wall -o nldft nldft.o
nldft.o: nldft.cpp
	c++ -Wall -c nldft.cpp
#maxwell_construction.o: maxwell_construction.cpp
#	c++ -Wall -c maxwell_construction.cpp

# make clean
clean:
	rm -f *.o nldft
