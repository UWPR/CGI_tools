CPP = g++
CPPFLAGS = -O2 -Wall -Wwrite-strings -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -fPIE

all: fragment

fragment: fragment_cgi.o cgi-util.o
	${CPP} fragment_cgi.o cgi-util.o -o fragment.cgi -lm

fragment_cgi.o: fragment_cgi.cpp AminoAcidMasses.h
	${CPP} ${CPPFLAGS} -c fragment_cgi.cpp -o fragment_cgi.o

cgi-util.o: cgi-util.cpp cgi-util.h
	${CPP} ${CPPFLAGS} -c cgi-util.cpp -o cgi-util.o

clean:
	rm -f *.o *.cgi
