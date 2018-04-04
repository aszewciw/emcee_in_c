
UTILS_DIR=..
# PARS_DIR=.
CC=mpicc
CFLAGS=-cc=icc -I$(UTILS_DIR) -Wall
CLINK= -limf

TARGET = run_mcmc
all: $(TARGET)

run_mcmc: test.o $(UTILS_DIR)/emcee.o
	$(CC) $(CFLAGS) -o $@ $^ $(CLINK)

clean:
	rm *.o
	rm $(UTILS_DIR)/*.o
	rm run_mcmc