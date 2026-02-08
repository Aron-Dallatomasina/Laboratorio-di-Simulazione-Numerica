CXX = g++
CXXFLAGS = -Wall -O3 -std=c++11 -I../include

main.exe: main.o random.o
	$(CXX) random.o main.o -o main.exe

main.o: main.cpp
	$(CXX) -c main.cpp -o main.o $(CXXFLAGS)

random.o: ../include/random.cpp ../include/random.h
	$(CXX) -c ../include/random.cpp -o random.o $(CXXFLAGS)

clean:
	rm -f *.o main.exe ../Output/seed.out
