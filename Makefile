CC=g++

CXXFLAGS=-Wall -pedantic -std=c++0x  -ggdb -DDEBUG

Optmizer=./optimizer
benchmark=./cec2013

#demoheader = $(wildcard ./*.h)
demosource = $(wildcard ./*.cpp)
demoobject = $(patsubst %.cpp, %.o, $(notdir $(demosource)))

cecheader=$(wildcard $(benchmark)/*.h)
cecsource=$(wildcard $(benchmark)/*.cpp)
cecobject=$(patsubst %.cpp, %.o, $(notdir $(cecsource)))

optheader=$(wildcard $(Optimizer)/*.h)
optsource=$(wldcard $(Optimizer)/*.cpp)
optobject=$(patsubst %.cpp, %.o, $(notdir $(optsource)))

OBJECTS=$(demoobject) $(cecobject) $(Optmizer)/CCJaDE.o

demo: $(OBJECTS)
	$(CC) $(CXXFLAGS) -o demo $(OBJECTS)

demo.o: $(demosource) $(cecheader) $(optheader) $(optsource)
	$(CC) $(CXXFLAGS) -c demo.cpp

SaNSDE.o: $(cecheader) $(Optmizer)/CCJaDE.h $(Optmizer)/CCJaDE.cpp
	$(CC) $(CXXFLAGS) -c $(optsource) 
	
Benchmarks.o:  $(benchmark)/Benchmarks.h $(benchmark)/Benchmarks.cpp
	$(CC) $(CXXFLAGS) -c $(benchmark)/Benchmarks.cpp

Header.o:  $(cecheader) $(cecsource)
	$(CC) $(CXXFLAGS) -c $(benchmark)/Header.cpp

F1.o: $(benchmark)/F1.h $(benchmark)/Benchmarks.h $(benchmark)/F1.cpp
	$(CC) $(CXXFLAGS) -c $(benchmark)/F1.cpp

F2.o: $(benchmark)/F2.h $(benchmark)/Benchmarks.h $(benchmark)/F2.cpp
	$(CC) $(CXXFLAGS) -c $(benchmark)/F2.cpp

F3.o: $(benchmark)/F3.h $(benchmark)/Benchmarks.h $(benchmark)/F3.cpp
	$(CC) $(CXXFLAGS) -c $(benchmark)/F3.cpp

F4.o: $(benchmark)/F4.h $(benchmark)/Benchmarks.h $(benchmark)/F4.cpp
	$(CC) $(CXXFLAGS) -c $(benchmark)/F4.cpp

F5.o: $(benchmark)/F5.h $(benchmark)/Benchmarks.h $(benchmark)/F5.cpp
	$(CC) $(CXXFLAGS) -c $(benchmark)/F5.cpp

F6.o: $(benchmark)/F6.h $(benchmark)/Benchmarks.h $(benchmark)/F6.cpp
	$(CC) $(CXXFLAGS) -c $(benchmark)/F6.cpp

F7.o: $(benchmark)/F7.h $(benchmark)/Benchmarks.h $(benchmark)/F7.cpp
	$(CC) $(CXXFLAGS) -c $(benchmark)/F7.cpp

F8.o: $(benchmark)/F8.h $(benchmark)/Benchmarks.h $(benchmark)/F8.cpp
	$(CC) $(CXXFLAGS) -c $(benchmark)/F8.cpp

F9.o: $(benchmark)/F9.h $(benchmark)/Benchmarks.h $(benchmark)/F9.cpp
	$(CC) $(CXXFLAGS) -c $(benchmark)/F9.cpp

F10.o: $(benchmark)/F10.h $(benchmark)/Benchmarks.h $(benchmark)/F10.cpp
	$(CC) $(CXXFLAGS) -c $(benchmark)/F10.cpp

F11.o: $(benchmark)/F11.h $(benchmark)/Benchmarks.h $(benchmark)/F11.cpp
	$(CC) $(CXXFLAGS) -c $(benchmark)/F11.cpp

F12.o: $(benchmark)/F12.h $(benchmark)/Benchmarks.h $(benchmark)/F12.cpp
	$(CC) $(CXXFLAGS) -c $(benchmark)/F12.cpp

F13.o: $(benchmark)/F13.h $(benchmark)/Benchmarks.h $(benchmark)/F13.cpp
	$(CC) $(CXXFLAGS) -c $(benchmark)/F13.cpp

F14.o: $(benchmark)/F14.h $(benchmark)/Benchmarks.h $(benchmark)/F14.cpp
	$(CC) $(CXXFLAGS) -c $(benchmark)/F14.cpp

F15.o: $(benchmark)/F15.h $(benchmark)/Benchmarks.h $(benchmark)/F15.cpp
	$(CC) $(CXXFLAGS) -c $(benchmark)/F15.cpp

.PHONY : clean
clean:
	rm -f demo $(OBJECTS)
