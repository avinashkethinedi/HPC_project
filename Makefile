CFLAGS = -O3 -lm
CC = gcc
MCC = mpicc

all: build

build: SpMV

SpMV: SpMV.o
	$(MCC) $(CFLAGS) $< -o $@

SpMV.o: SpMV.c
	$(MCC) $(CFLAGS) $< -c -o $@

clean:
	rm -rf SpMV core* *.o
