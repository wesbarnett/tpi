
PKGNAME = tpi
PREFIX ?= /usr/local

.PHONY: all install clean

all:
	mkdir -p bin
	g++ src/main.cpp src/Atomtype.cpp src/Ini.cpp -O2 -o bin/tpi -I./include -fopenmp -Wall -march=native -lgmxcpp

install:
	install -Dm755 bin/* -t $(DESTDIR)$(PREFIX)/bin/

clean:
	rm -r bin
