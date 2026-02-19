CXX = mpic++
CXXFLAGS = -Wall -O3 -std=c++17 -I../Include

main.exe: main.o random.o genetic.o
	$(CXX) main.o random.o genetic.o -o main.exe

main.o: main.cpp
	$(CXX) -c main.cpp -o main.o $(CXXFLAGS)

random.o: ../Include/random.cpp ../Include/random.h
	$(CXX) -c ../Include/random.cpp -o random.o $(CXXFLAGS)

genetic.o: ../Include/genetic.cpp ../Include/genetic.h
	$(CXX) -c ../Include/genetic.cpp -o genetic.o $(CXXFLAGS)

clean:
	rm -f *.o main.exe




