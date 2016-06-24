all:
	g++ main.cpp Atomtype.cpp -O2 /home/wes/libgmxcpp/lib/libgmxcpp.so -o run -I. -I/home/wes/libgmxcpp/include -fopenmp -Wall -march=native
