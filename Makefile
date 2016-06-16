all:
	g++ main.cpp -O2 Atomtype.cpp /home/wes/libgmxcpp/lib/libgmxcpp.so -o run -I. -I/home/wes/libgmxcpp/include -fopenmp -Wall -march=native
