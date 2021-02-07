#!/bash/bin

CXX = g++
CXXFLAGS = -std=c++11 -O2 -lm  -lgsl -lgslcblas -fopenmp -lfftw3
objdir=obj
srcs=functions.cpp main.cpp 
objects = functions.o main.o 

prog: main.o functions.o  WB.o
	$(CXX)  -o $@  $^ $(CXXFLAGS)

main.o : main.cpp functions.h WB.h 
	$(CXX) -c main.cpp  $(CXXFLAGS)

functions.o: functions.cpp functions.h  WB.h
	$(CXX) -c functions.cpp $(CXXFLAGS)

WB.o: WB.cpp functions.h  WB.h
	$(CXX) -c WB.cpp $(CXXFLAGS)
clean:
	rm -rf  *.pyc a.out ./../data/text/*.txt ./../data/*.png 
	rm -f *.png *.out ../data/*.bin ../data/*.pdf
	rm -f ../data/fig/* ../data/npz/* 
eradicate:
	rm -rf *.o  *.pyc a.out ../data/text/*.txt  ../data/*.png 
	rm -f *.png prog *.out ../data/*.dat ../data/*.bin
	rm -f ../data/*.pdf
	rm -f ../data/fig/* ../data/npz/* 

.PHONY: clean, eradicate
