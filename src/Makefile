PS_PATH = -I../fdps

#CC = g++
#CC = mpic++
CC = /opt/local/bin/mpic++-mpich-mp
CFLAGS = -std=c++11
CFLAGS += -O3
CFLAGS += -Wall
CFLAGS += -ffast-math
CFLAGS += -funroll-loops

#CFLAGS += -g -O0

#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

CFLAGS += -DNOSE_HOOVER
#CFLAGS += -DDISSIPATIVE_RANDOM

CFLAGS += -DSMALL_SOLVENT

CFLAGS += -DNANOSLIT
#CFLAGS += -DNANOTUBE


SRC = janus.cpp
PROGRAM = $(SRC:%.cpp=%.out)

$(PROGRAM):$(SRC)
	$(CC) $(MULEXP) $(PS_PATH) $(CFLAGS) -o $@ $<

clean:
	rm $(PROGRAM)


