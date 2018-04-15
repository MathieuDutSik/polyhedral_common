example: example.o
	g++ -std=c++11 example.o -o example $(MPFR_LINK) -lgmp

example.o: example.cpp
	g++ -std=c++11 -I$(MPFR_INCDIR) -c example.cpp

clean:
	rm -f *.o 

