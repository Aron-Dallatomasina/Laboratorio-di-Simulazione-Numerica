CXX = g++
CXXFLAGS = -Wall -O3 -std=c++11 -I../Include

Monte_Carlo_ISING_1D.exe: Monte_Carlo_ISING_1D.o random.o
	$(CXX) random.o Monte_Carlo_ISING_1D.o -o Monte_Carlo_ISING_1D.exe

Monte_Carlo_ISING_1D.o: Monte_Carlo_ISING_1D.cpp Monte_Carlo_ISING_1D.h
	$(CXX) -c Monte_Carlo_ISING_1D.cpp -o Monte_Carlo_ISING_1D.o $(CXXFLAGS)

random.o: ../Include/random.cpp ../Include/random.h
	$(CXX) -c ../Include/random.cpp -o random.o $(CXXFLAGS)

clean:
	rm -f *.o Monte_Carlo_ISING_1D.exe