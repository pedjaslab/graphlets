CC=g++
DEBUG=-g
FLAGS=-Wall -O2
LIBS=-lm
OBJS=utils.o string.o


all: run_kernel



utils.o: utils.h utils.cc config.h
	$(CC) $(DEBUG) $(FLAGS) -c utils.cc

string.o: string.h string.cc
	$(CC) $(DEBUG) $(FLAGS) -c string.cc

simgraph.o: simgraph.h simgraph.cc string.h
	$(CC) $(DEBUG) $(FLAGS) -c simgraph.cc




digkernel.o: digkernel.h digkernel.cc config.h
	$(CC) $(DEBUG) $(FLAGS) -c digkernel.cc

run_kernel: run_kernel.cc digkernel.o simgraph.o $(OBJS) 
	$(CC) $(DEBUG) $(FLAGS) $(LIBS) run_kernel.cc digkernel.o simgraph.o $(OBJS) -o $@



clean:
	rm -rf run_kernel *.o *.*~ *~ core *.dSYM

