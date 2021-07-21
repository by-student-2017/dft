# make
#nldft: qsdft.o maxwell_construction.o
#	c++ -Wall -o qsdft.o maxwell_construction.o
nldft: nldft.o
	c++ -Wall -o qsdft qsdft.o
nldft.o: nldft.cpp
	c++ -Wall -c qsdft.cpp
#maxwell_construction.o: maxwell_construction.cpp
#	c++ -Wall -c maxwell_construction.cpp

# make clean
clean:
	rm -f *.o nldft
