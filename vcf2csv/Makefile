CPP=g++
BOOST_ROOT=/opt/boost-1.54.0
CPPFLAGS=--std=c++11 -ggdb -Wall -Wno-unused-local-typedefs -O2 -I$(BOOST_ROOT)/include
LDFLAGS=-L$(BOOST_ROOT)/lib64
BOOST_LIBS=-Wl,-R$(BOOST_ROOT)/lib64 -lboost_program_options 

all: vcf2csv

vcf2csv: vcf2csv.cpp
	$(CPP) $(CPPFLAGS) $(LDFLAGS) -o $@ $<

clean:
	rm vcf2csv

