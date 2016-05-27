all:
	g++ main.cpp Atomtype.cpp -o run -I. -lgmxcpp -fopenmp -Wall
