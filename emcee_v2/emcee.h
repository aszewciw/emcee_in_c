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

/* position of a single walker */
typedef struct walker_pos{
    int accept;         /* were the current pars part of an accepted step or not (1 or 0) */
    double lnprob;      /* ln of probability for current pars */
    double pars[NPARS]; /* parameters */
} walker_pos;

/* an ensemble of walkers */
typedef struct ensemble {
    size_t nwalkers;    /* number of walkers in ensemble */
    size_t npars;       /* number of pars */
    walker_pos *walker; /* positions of walkers in ensemble */
} ensemble;

/* stored data for the full chain */
typedef struct chain {
    size_t nsteps;          /* number of steps in chain */
    ensemble * ensemble_A;  /* first ensemble of walkers */
    ensemble * ensemble_B;  /* second ensemble of walkers */
} chain;

/* allocate space for walkers, ensembles, chains */
walker_pos* allocate_walkers(size_t nwalkers);
ensemble* allocate_ensemble(size_t nwalkers, size_t npars);
chain* allocate_chain(size_t nsteps, size_t nwalkers, size_t npars);

/* free space for walkers, ensembles, chains */
void free_walkers(walker_pos *w);
void free_ensemble(ensemble *e);
void free_chain(chain *c);


walker_pos* make_guess(double *centers, double *widths, size_t nwalkers, size_t npars);

/* generate random number between 0 and 1 */
double rand_0to1();

/* generate normal random numbers */
double normal_rand();

/* decide to accept or reject proposed parameters based on new and old probabilities */
int walker_accept(double lnprob_old,double lnprob_new,size_t npars,double z);

/* write every step in mcmc */
int write_chain(const struct chain *c, const char *fname);

/* write single step in chain */
int write_step(const char *fname, const struct ensemble *ensemble_A, const struct ensemble *ensemble_B);

/* write header info of chain */
int write_header(const char *fname, size_t nsteps, size_t nwalkers, size_t npars);

/* use positions of walkers in the complementary ensemble to move a single walker */
void step_walkers(walker_pos *walkers, ensemble *comp_walkers, size_t nwalkers,
                  double a, double (*lnprob)(const double *, size_t, const void *),
                  const void *userdata);

/* choose a random integer between 0 and n */
size_t rand_walker(size_t n);

/* get a random z for walker step from distribution g(z) */
double rand_gofz(double a);

/* the main function to run a chain */
void run_chain(int *argc, char ***argv, walker_pos *start_pos, double a,
               double (*lnprob)(const double *, size_t, const void *),
               const void *userdata, const char *fname, int nburn_in);

#endif //#ifndef _EMCEE_HEADER_GUARD
