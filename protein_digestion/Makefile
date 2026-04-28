CPP = g++
CPPFLAGS = -O2 -Wall -Wwrite-strings -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -fPIE

all: digest

digest: digest_cgi.o cgi-util.o
	${CPP} digest_cgi.o cgi-util.o -o digest.cgi -lm

digest_cgi.o: digest_cgi.c
	${CPP} ${CPPFLAGS} -c digest_cgi.c -o digest_cgi.o

cgi-util.o: cgi-util.cpp cgi-util.h
	${CPP} ${CPPFLAGS} -c cgi-util.cpp -o cgi-util.o

clean:
	rm -f *.o *.cgi
