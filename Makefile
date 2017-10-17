

CC=gcc
CFLAGS=-std=gnu99 -Wall -Werror -O2
LFLAGS=-lm

#################################### compile ###################################
all: mca-const mca-rand

mca-const: test-mca-const.o mca.o rand.o
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

mca-rand: test-mca-rand.o mca.o rand.o
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

clean:
	rm -f *.o
	rm -f mca-const
	rm -f mca-rand