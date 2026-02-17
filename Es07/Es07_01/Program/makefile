CXX = g++
CXXFLAGS = -Wall -O3 -std=c++11 -I../Include

NVE_NVT.exe: MD_MC.o random.o
	$(CXX) random.o MD_MC.o -o NVE_NVT.exe

MD_MC.o: MD_MC.cpp MD_MC.h
	$(CXX) -c MD_MC.cpp -o MD_MC.o $(CXXFLAGS)

random.o: ../Include/random.cpp ../Include/random.h
	$(CXX) -c ../Include/random.cpp -o random.o $(CXXFLAGS)

clean:
	rm -f *.o NVE_NVT.exe

