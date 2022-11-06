CFLAGS = -O3 -lm
CC = gcc
MCC = mpic++

all: build

build: SpMV

SpMV: SpMV.o
	$(MCC) $(CFLAGS) $< -o $@

SpMV.o: SpMV.cpp
	$(MCC) $(CFLAGS) $< -c -o $@

clean:
	rm -rf SpMV core* *.o
