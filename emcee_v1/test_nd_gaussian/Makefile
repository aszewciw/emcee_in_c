
UTILS_DIR=..
PARS_DIR=.
CC=mpicc
CFLAGS=-cc=icc -I$(UTILS_DIR) -I$(PARS_DIR) -Wall

TARGET = run_mcmc
all: $(TARGET)

run_mcmc: test.o $(UTILS_DIR)/emcee.o
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm *.o
	rm $(UTILS_DIR)/*.o
	rm run_mcmc