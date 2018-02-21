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

#define WORKTAG 1
#define DIETAG 2

/* position of a single walker */
typedef struct walker_pos{
    int accept;         /* were the current pars part of an accepted step or not (1 or 0) */
    double lnprob;      /* ln of probability for current pars */
    double pars[NPARS]; /* parameters */
} walker_pos;

// /* an ensemble of walkers */
// typedef struct ensemble {
//     size_t nwalkers;    /* number of walkers in ensemble */
//     size_t npars;       /* number of pars */
//     walker_pos *walker; /* positions of walkers in ensemble */
// } ensemble;

//  stored data for the full chain
// typedef struct chain {
//     size_t nsteps;          /* number of steps in chain */
//     ensemble * ensemble_A;  /* first ensemble of walkers */
//     ensemble * ensemble_B;  /* second ensemble of walkers */
// } chain;

/* allocate space for walkers, ensembles, chains */
walker_pos* allocate_walkers(size_t nwalkers);
// ensemble* allocate_ensemble(size_t nwalkers, size_t npars);
// chain* allocate_chain(size_t nsteps, size_t nwalkers, size_t npars);

/* free space for walkers, ensembles, chains */
void free_walkers(walker_pos *w);
// void free_ensemble(ensemble *e);
// void free_chain(chain *c);


walker_pos* make_guess(double *centers, double *widths, size_t nwalkers, size_t npars);

/* generate random number between 0 and 1 */
double rand_0to1();

/* generate normal random numbers */
double normal_rand();

/* decide to accept or reject proposed parameters based on new and old probabilities */
int walker_accept(double lnprob_old, double lnprob_new, size_t npars, double z);

/* write every step in mcmc */
// int write_chain(const struct chain *c, const char *fname);

/* write single step in chain */
int write_step(const char *fname, const struct walker_pos *ensemble_A,
               const struct walker_pos *ensemble_B, size_t nwalkers, size_t istep);

/* write header info of chain */
int write_header(const char *fname, size_t nsteps, size_t nwalkers, size_t npars);


void create_trials(walker_pos *trial, const struct walker_pos *walkers,
                   const struct walker_pos *comp_walkers, double *z_array,
                   size_t nwalkers, double a);

void update_positions(walker_pos *walkers, const struct walker_pos *trial,
                      double *z_array, size_t nwalkers);

void manager(walker_pos *start_pos, double a, const char *fname, int nburn);

void worker(const void *userdata, double (*lnprob)(const double *, size_t, const void *));

/* choose a random integer between 0 and n */
size_t rand_walker(size_t n);

/* get a random z for walker step from distribution g(z) */
double rand_gofz(double a);

// void run_chain_standard(int *argc, char ***argv, walker_pos *start_pos, double a,
//                         double (*lnprob)(const double *, size_t, const void *),
//                         const void *userdata, const char *fname, int nburn);

void run_chain_loadbalancing(int *argc, char ***argv, walker_pos *start_pos, double a,
                             double (*lnprob)(const double *, size_t, const void *),
                             const void *userdata, const char *fname, int nburn);

/* the main function to run a chain */
void run_chain(int *argc, char ***argv, walker_pos *start_pos, double a,
               double (*lnprob)(const double *, size_t, const void *),
               const void *userdata, const char *fname, int nburn,
               int load_balancing);

#endif //#ifndef _EMCEE_HEADER_GUARD
