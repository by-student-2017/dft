# make
#qsdft: qsdft_fmt_slit.o maxwell_construction.o
#	c++ -O2 -Wall -o qsdft_fmt_slit.o maxwell_construction.o
qsdft: qsdft.o
	c++ -O2 -Wall -o qsdft_fmt_slit qsdft_fmt_slit.o
qsdft.o: qsdft_fmt_slit.cpp
	c++ -O2 -Wall -c qsdft_fmt_slit.cpp
#maxwell_construction.o: maxwell_construction.cpp
#	c++ -O2 -Wall -c maxwell_construction.cpp

# make clean
clean:
	rm -f *.o nldft
