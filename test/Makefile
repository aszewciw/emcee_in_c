

INC_DIR=../
CC=gcc
CFLAGS=-std=gnu99 -Wall -Werror -O2 -I $(INC_DIR)
LFLAGS=-lm

#################################### compile ###################################
all: test-mca-const test-mca-rand

mca-const: test-mca-const.o $(INC_DIR)/mca.o $(INC_DIR)/rand.o
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

mca-rand: test-mca-rand.o $(INC_DIR)/mca.o $(INC_DIR)/rand.o
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

clean:
	rm -f *.o
	rm -f ../*.o
	rm -f mca-const
	rm -f mca-rand