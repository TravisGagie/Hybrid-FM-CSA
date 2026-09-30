CXX ?= g++
CXXFLAGS ?= -O3 -march=native -DNDEBUG -std=c++17
hybrid: hybrid.cpp runrank.hpp efrank.hpp huffwt.hpp libsais/src/libsais.c
	gcc -O3 -march=native -DNDEBUG -c libsais/src/libsais.c -Ilibsais/include -o libsais.o
	$(CXX) $(CXXFLAGS) -Ilibsais/include hybrid.cpp libsais.o -o hybrid
clean:
	rm -f hybrid libsais.o
