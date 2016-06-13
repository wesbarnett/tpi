all:
	g++ main.cpp -Ofast Atomtype.cpp -o run -I. -lgmxcpp -fopenmp -Wall
