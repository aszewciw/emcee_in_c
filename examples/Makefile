
EMCEE_DIR=..
CC=mpicc
CFLAGS=-cc=icc -I$(UTILS_DIR) -Wall
CLINK= -limf

TARGET = run_mcmc
all: $(TARGET)

run_mcmc: nd_gaussian.o $(EMCEE_DIR)/emcee.o
	$(CC) $(CFLAGS) -o $@ $^ $(CLINK)

clean:
	rm *.o
	rm $(EMCEE_DIR)/*.o
	rm run_mcmc