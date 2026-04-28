# Makefile for motif_search CGI + test suite
#
# Targets:
#   make              — build motif_search.cgi
#   make test         — build everything and run all tests
#   make clean        — remove compiled binaries

CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

.PHONY: all test clean

all: motif_search.cgi

motif_search.cgi: motif_search.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

motif_search_tests: motif_search_tests.cpp motif_search.cgi
	$(CXX) $(CXXFLAGS) -o $@ motif_search_tests.cpp

test: motif_search_tests
	./motif_search_tests

clean:
	rm -f motif_search.cgi motif_search_tests /tmp/cgi_test_in.bin /tmp/cgi_test_out.txt
