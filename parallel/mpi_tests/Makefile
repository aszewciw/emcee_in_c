CC=mpicc
CFLAGS=-Wall

step_comm : step_comm.c
	$(CC) -o $@ $< ${CFLAGS}

.PHONY: clean

clean:
	rm -f step_comm