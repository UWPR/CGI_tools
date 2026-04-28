CPP = g++
CPPFLAGS = -O2 -Wall -Wwrite-strings -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -fPIE

all: element

element: element_cgi.o cgi-util.o
	${CPP} element_cgi.o cgi-util.o -o element.cgi -lm

element_cgi.o: element_cgi.c AminoAcidMasses.h
	${CPP} ${CPPFLAGS} -c element_cgi.c -o element_cgi.o

cgi-util.o: cgi-util.cpp cgi-util.h
	${CPP} ${CPPFLAGS} -c cgi-util.cpp -o cgi-util.o

clean:
	rm -f *.o *.cgi
