all:
	g++ main.cpp -O2 /home/wes/libgmxcpp/lib/libgmxcpp.so -o run -I/home/wes/libgmxcpp/include -fopenmp -Wall -march=native
