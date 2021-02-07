#!/bash/bin

CXX = g++
CXXFLAGS = -std=c++11 -O2 -lm  -lgsl -lgslcblas -fopenmp -lfftw3
objdir=obj

objects = lib.o main.o 

prog: main.o lib.o  Traub.o
	$(CXX)  -o $@  $^ $(CXXFLAGS)

main.o : main.cpp lib.h Traub.h
	$(CXX) -c main.cpp  $(CXXFLAGS)

lib.o: lib.cpp lib.h  Traub.h
	$(CXX) -c lib.cpp $(CXXFLAGS)

Traub.o: Traub.cpp lib.h  Traub.h
	$(CXX) -c Traub.cpp $(CXXFLAGS)
clean:
	rm -rf  *.pyc 
	rm -rf ./../data/text/*
	rm -rf ../data/fig/*
	rm -rf ../data/npz/* 

eradicate: clean
	rm -rf *.o   
	rm -f prog *.out
	

.PHONY: clean, eradicate
