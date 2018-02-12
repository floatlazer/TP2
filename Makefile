# Makefile TD n°2
#CXX = g++ -std=c++11
CXX = mpic++ -std=c++11
CXXFLAGS = -O2 -Wall -pedantic -march=native
#CXXFLAGS = -g -Wall -pedantic -march=native
LIBS = -lm

all:    Mandelbrot.exe  Mandelbrot_SepLine.exe Mandelbrot_MasterSlave.exe matvec.exe


Mandelbrot.exe: Mandelbrot.cpp lodepng/lodepng.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

Mandelbrot_SepLine.exe: Mandelbrot_SepLine.cpp lodepng/lodepng.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

Mandelbrot_MasterSlave.exe: Mandelbrot_MasterSlave.cpp lodepng/lodepng.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

matvec.exe:	matvec.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

cleanall:
	@rm -rf *.o *~ *.exe
