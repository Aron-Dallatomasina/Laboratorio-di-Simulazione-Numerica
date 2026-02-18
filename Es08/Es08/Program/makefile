CXX = g++
CXXFLAGS = -Wall -O3 -std=c++11 -I../Include

main.exe: main.o random.o
	$(CXX) random.o main.o -o main.exe

main.o: main.cpp
	$(CXX) -c main.cpp -o main.o $(CXXFLAGS)

random.o: ../Include/random.cpp ../Include/random.h
	$(CXX) -c ../Include/random.cpp -o random.o $(CXXFLAGS)

clean:
	rm -f *.o main.exe

