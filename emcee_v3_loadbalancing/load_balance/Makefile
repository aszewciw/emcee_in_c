
CC=mpicc
CFLAGS=-cc=icc -Wall
CLINK= -limf

TARGET = load_bal
all: $(TARGET)

load_bal: test.o
	$(CC) $(CFLAGS) -o $@ $^ $(CLINK)

clean:
	rm *.o
	rm load_bal