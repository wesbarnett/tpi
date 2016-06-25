all:
	g++ src/main.cpp src/Atomtype.cpp src/Ini.cpp -O2 /home/wes/libgmxcpp/lib/libgmxcpp.so -o run -I./include -I/home/wes/libgmxcpp/include -fopenmp -Wall -march=native
