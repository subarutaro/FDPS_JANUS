PS_PATH = -I../fdps

CC = g++
#CC = mpicxx
CFLAGS = -std=c++11
CFLAGS += -O3
CFLAGS += -Wall
CFLAGS += -ffast-math
CFLAGS += -funroll-loops
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

#CFLAGS += -DNOSE_HOOVER
CFLAGS += -DDSSIPATIVE_RANDOM

#CFLAGS += -DNANOSLIT
CFLAGS += -DNANOTUBE



SRC = janus.cpp
PROGRAM = $(SRC:%.cpp=%.out)

$(PROGRAM):$(SRC)
	$(CC) $(MULEXP) $(PS_PATH) $(CFLAGS) -o $@ $<

a.out: test.cpp
	$(CC) $(MULEXP) $(PS_PATH) $(CFLAGS) test.cpp

clean:
	rm $(PROGRAM)


