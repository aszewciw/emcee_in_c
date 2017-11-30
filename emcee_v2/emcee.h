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

typedef struct walker_pos{
    int accept;
    double lnprob;
    double pars[NPARS];
} walker_pos;

typedef struct ensemble {
    size_t nwalkers;
    size_t npars;
    walker_pos *walker;
} ensemble;

typedef struct chain {
    size_t nsteps;
    ensemble * ensemble_A;
    ensemble * ensemble_B;
} chain;

walker_pos* allocate_walkers(size_t nwalkers);
ensemble* allocate_ensemble(size_t nwalkers, size_t npars);
chain* allocate_chain(size_t nsteps, size_t nwalkers, size_t npars);

void free_walkers(walker_pos *w);
void free_ensemble(ensemble *e);
void free_chain(chain *c);

walker_pos* make_guess(double *centers, double *widths, size_t nwalkers, size_t npars);
double rand_0to1();
double normal_rand();

int walker_accept(double lnprob_old,double lnprob_new,size_t npars,double z);
int write_chain(const struct chain *c, const char *fname);
void step_walkers(walker_pos *walkers, ensemble *comp_walkers, size_t nwalkers,
                  double a, double (*lnprob)(const double *, size_t, const void *),
                  const void *userdata);

size_t rand_walker(size_t n);
double rand_gofz(double a);
void run_chain(int *argc, char ***argv, walker_pos *start_pos, double a,
               double (*lnprob)(const double *, size_t, const void *),
               const void *userdata, const char *fname);

#endif //#ifndef _EMCEE_HEADER_GUARD
