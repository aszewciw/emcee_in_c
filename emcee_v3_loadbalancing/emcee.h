#ifndef _EMCEE_HEADER_GUARD
#define _EMCEE_HEADER_GUARD

// #include "pars.h"
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define WORKTAG (1)
#define DIETAG (2)

/* position of a single walker */
typedef struct walker_pos{
    int rank;           /* process rank that handled the computation */
    int accept;         /* were the current pars part of an accepted step or not (1 or 0) */
    double lnprob;      /* ln of probability for current pars */
    // double pars[NPARS]; /* parameters */
    double *pars; /* parameters */
} walker_pos;

/* get number of lines in a file */
size_t getNlines(const char *fname,const char comment);

/* allocate space for walkers */
walker_pos* allocate_walkers(int nwalkers, int npars);

/* free space for walkers */
void free_walkers(int nwalkers, walker_pos *w);

/* creates a "ball" of initial guesses (walker positions) */
walker_pos *make_guess(int nwalkers, int npars, double *centers, double *widths);

/* decide to accept or reject proposed parameters based on new and old probabilities */
int walker_accept(int npars, double lnprob_old, double lnprob_new, double z)

/* generate random number between 0 and 1 */
double rand_0to1(void);

/* generate normal random numbers */
double normal_rand(void);

/* choose a random integer between 0 and max */
unsigned long rand_walker(unsigned long max);

/* get a random z for walker step from distribution g(z) */
double rand_gofz(double a);

/* write header info of chain */
int write_header(int nsteps, int nwalkers, int npars, const char *fname)

/* write single step in chain */
int write_step(int npars, int nwalkers, int istep, const struct walker_pos *ensemble_A,
               const struct walker_pos *ensemble_B, const char *fname)

/* creates an ensemble of "trial" walker positions */
void create_trials(int nwalkers, int npars, double a, double *z_array,
                   walker_pos *trial, const struct walker_pos *walkers,
                   const struct walker_pos *comp_walkers);

/* accept or reject trial positions */
void update_positions(int nwalkers, int npars, double *z_array, walker_pos *walkers,
                      const struct walker_pos *trial);

/*
  manager (process 0) is responsible for:
    (1) storing and updating the state of the ensemble
    (2) choosing new trial positions
    (3) distibuting trial positions among workers (when they are ready)
    (4) accepting/rejecting trial positions
    (5) writing the ensemble state to a file
 */
void manager(int nwalkers, int nsteps, int npars, int nburn, int resume, double a,
             walker_pos *start_pos, const char *fname);

/*
  worker (processes 1 - nprocs) is responsible for:
    (1) accepting parameters from manager
    (2) computing lnprob for walker positions
 */
void worker(int npars, const void *userdata, double (*lnprob)(const double *, int, const void *));

/* run chain with loadbalancing enabled */
void run_chain(int *argc, char ***argv, int nwalkers, int nsteps, int npars,
               int nburn, int resume, double a, walker_pos *start_pos,
               double (*lnprob)(const double *, int, const void *),
               const void *userdata, const char *fname);

#endif //#ifndef _EMCEE_HEADER_GUARD
