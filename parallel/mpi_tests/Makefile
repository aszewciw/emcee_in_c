CC=mpicc
CFLAGS=-Wall

step_comm_alt : step_comm_alt.c
	$(CC) -o $@ $< ${CFLAGS}

.PHONY: clean

clean:
	rm -f step_comm_alt