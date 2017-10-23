CC=mpicc
CFLAGS=-Wall

struct_comm : struct_comm.c
	$(CC) -o $@ $< ${CFLAGS}

.PHONY: clean

clean:
	rm -f student
	rm -f struct_comm