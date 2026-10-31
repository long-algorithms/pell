CXX ?= g++
CXXFLAGS ?= -std=c++17 -O3 -DNDEBUG -march=native -Wall -Wextra
LDLIBS ?= -lgmpxx -lgmp

all: conway_fast

conway_fast: conway_fast.cpp primes_10000.hpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDLIBS)

test: conway_fast check-square-factor
	python3 test_cli.py

clean:
	rm -f conway_fast check-square-factor

.PHONY: all test clean

check-square-factor: test_square_factor.cpp conway_fast.cpp primes_10000.hpp
	$(CXX) $(CXXFLAGS) test_square_factor.cpp -o check-square-factor $(LDLIBS)
	./check-square-factor
