CXX = g++ -O3 -Wall -Wwrite-strings
CXXFLAGS = -c -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64


# Link command:
all: calcisotopes

calcisotopes: calcisotopes_cgi.o  isotopes.o cgi-util.o
	${CXX} calcisotopes_cgi.o isotopes.o cgi-util.o -o calcisotopes.cgi -lm

calcisotopes_cgi.o: calcisotopes_cgi.cpp
	${CXX} ${CXXFLAGS} calcisotopes_cgi.cpp -o calcisotopes_cgi.o

isotopes.o: isotopes.cpp
	${CXX} ${CXXFLAGS} isotopes.cpp -o isotopes.o

cgi-util.o: cgi-util.cpp
	${CXX} ${CXXFLAGS} cgi-util.cpp -o cgi-util.o

clean:
	rm -f *.o *.cgi
