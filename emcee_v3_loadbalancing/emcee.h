#ifndef _EMCEE_HEADER_GUARD
#define _EMCEE_HEADER_GUARD

#include "pars.h"
#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define WORKTAG (1)
#define DIETAG (2)

/* position of a single walker */
typedef struct walker_pos{
    int rank;           /* process rank that handled the computation */
    int accept;         /* were the current pars part of an accepted step or not (1 or 0) */
    double lnprob;      /* ln of probability for current pars */
    double pars[NPARS]; /* parameters */
} walker_pos;

/* get number of lines in a file */
size_t getNlines(const char *fname,const char comment);

/* allocate space for walkers */
walker_pos* allocate_walkers(size_t nwalkers);

/* free space for walkers */
void free_walkers(walker_pos *w);

/* creates a "ball" of initial guesses (walker positions) */
walker_pos* make_guess(double *centers, double *widths, size_t nwalkers, size_t npars);

/* generate random number between 0 and 1 */
double rand_0to1();

/* generate normal random numbers */
double normal_rand();

/* decide to accept or reject proposed parameters based on new and old probabilities */
int walker_accept(double lnprob_old, double lnprob_new, size_t npars, double z);

/* write single step in chain */
int write_step(const char *fname, const struct walker_pos *ensemble_A,
               const struct walker_pos *ensemble_B, size_t nwalkers, size_t istep);

/* write header info of chain */
int write_header(const char *fname, size_t nsteps, size_t nwalkers, size_t npars);

/* creates an ensemble of "trial" walker positions */
void create_trials(walker_pos *trial, const struct walker_pos *walkers,
                   const struct walker_pos *comp_walkers, double *z_array,
                   size_t nwalkers, double a);

/* accept or reject trial positions */
void update_positions(walker_pos *walkers, const struct walker_pos *trial,
                      double *z_array, size_t nwalkers);

/*
  manager (process 0) is responsible for:
    (1) storing and updating the state of the ensemble
    (2) choosing new trial positions
    (3) distibuting trial positions among workers (when they are ready)
    (4) accepting/rejecting trial positions
    (5) writing the ensemble state to a file
 */
void manager(walker_pos *start_pos, double a, const char *fname, int nburn, int resume);

/*
  worker (processes 1 - nprocs) is responsible for:
    (1) accepting parameters from manager
    (2) computing lnprob for walker positions
 */
void worker(const void *userdata, double (*lnprob)(const double *, size_t, const void *));

/* choose a random integer between 0 and n */
size_t rand_walker(size_t n);

/* get a random z for walker step from distribution g(z) */
double rand_gofz(double a);

/* run chain with loadbalancing enabled */
void run_chain_loadbalancing(int *argc, char ***argv, walker_pos *start_pos, double a,
                             double (*lnprob)(const double *, size_t, const void *),
                             const void *userdata, const char *fname, int nburn, int resume);

/* the main function to run a chain */
void run_chain(int *argc, char ***argv, walker_pos *start_pos, double a,
               double (*lnprob)(const double *, size_t, const void *),
               const void *userdata, const char *fname, int nburn, int resume,
               int load_balancing);

#endif //#ifndef _EMCEE_HEADER_GUARD
