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
    int nwalkers;
    int npars;
    walker_pos *walker;
} ensemble;

typedef struct chain {
    int nsteps;
    ensemble * ball_1;
    ensemble * ball_2;
} chain;

walker_pos* allocate_walkers(int nwalkers);
ensemble* allocate_ensemble(int nwalkers, int npars);
chain* allocate_chain(int nsteps, int nwalkers, int npars);

void free_walkers(walker_pos *w);
void free_ensemble(ensemble *e);
void free_chain(chain *c);

struct walker_pos *make_guess(double *centers, double *widths, int nwalkers, int npars);
double rand_0to1();
double mca_randn();

int walker_accept(double lnprob_old,double lnprob_new,int npars,double z);
int write_chain(const struct chain *c, const char *fname);
void step_walkers(walker_pos *walkers, ensemble *comp_walkers, int nwalkers,
                  double a, double (*lnprob)(const double *, int, const void *),
                  const void *userdata);

int rand_walker(int n);
double rand_gofz(double a);

#endif //#ifndef _EMCEE_HEADER_GUARD
