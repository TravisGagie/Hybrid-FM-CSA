CXX ?= g++
CXXFLAGS ?= -O3 -march=native -DNDEBUG -std=c++17
psie: psie.cpp runrank.hpp efrank.hpp libsais/src/libsais.c
	gcc -O3 -march=native -DNDEBUG -c libsais/src/libsais.c -Ilibsais/include -o libsais.o
	$(CXX) $(CXXFLAGS) -Ilibsais/include psie.cpp libsais.o -o psie
clean:
	rm -f psie libsais.o
