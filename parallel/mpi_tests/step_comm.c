#include <stdio.h>
#include <stdlib.h>
#include <stddef.h> // offsetof()
#include "mpi.h"


// this is the custom structure that will be communicated in MPI
typedef struct walker_pos {
  int accept;
  double lnprob;
  double *pars;
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

// walker_pos* allocate_walker(){

// }

// ensemble* allocate_ensemble(){

// }

chain* allocate_chain(int nsteps, int nwalkers, int npars){

  int nwalkers_over_two=nwalkers/2;
  struct chain *self=calloc(1,sizeof(chain));
  if (self==NULL) {
      fprintf(stderr,"Could not allocate struct chain\n");
      exit(EXIT_FAILURE);
  }

  fprintf(stderr, "Address of chain nsteps: %p\n", (void*)&self->nsteps);
  fprintf(stderr, "Address of chain ball_1: %p\n", (void*)&self->ball_1);
  fprintf(stderr, "Address of chain ball_2: %p\n", (void*)&self->ball_2);


  /* allocate space for nsteps ensembles */
  self->ball_1=calloc(nsteps,sizeof(ensemble));
  if (self->ball_1==NULL) {
      fprintf(stderr,"Could not allocate struct ensemble\n");
      exit(EXIT_FAILURE);
  }
  self->ball_2=calloc(nsteps,sizeof(ensemble));
  if (self->ball_2==NULL) {
      fprintf(stderr,"Could not allocate struct ensemble\n");
      exit(EXIT_FAILURE);
  }

  for(int i=0; i<nsteps; i++){
    /*
    I realize its somewhat dumb to not make nwalkers part of the chain struct,
    but I'd prefer it to be a property of an ensemble, even though in this case
    I'm copying the same number nsteps times.

    Same can be said for npars. I could have also made the extra annoying decision
    to place it within the struct walker_pos. But this is the custom struct I
    set up for an MPI Allgatherv. It seems quite dumb to be gathering a value
    that never changes...
    */
    self->ball_1[i].nwalkers=nwalkers_over_two;
    self->ball_1[i].npars=npars;
    self->ball_1[i].walker=calloc(nwalkers_over_two,sizeof(walker_pos));
    fprintf(stderr, "Address of chain ball_1 walker_pos %d: %p\n", i,(void*)&self->ball_1[i].walker);
    if (self->ball_1[i].walker==NULL) {
        fprintf(stderr,"Could not allocate struct walker_pos\n");
        exit(EXIT_FAILURE);
    }
    self->ball_2[i].nwalkers=nwalkers_over_two;
    self->ball_2[i].npars=npars;
    self->ball_2[i].walker=calloc(nwalkers_over_two,sizeof(walker_pos));
    if (self->ball_2[i].walker==NULL) {
        fprintf(stderr,"Could not allocate struct walker_pos\n");
        exit(EXIT_FAILURE);
    }
    for(int j=0; j<nwalkers_over_two; j++){
      self->ball_1[i].walker[j].pars=calloc(npars,sizeof(double));
      if (self->ball_1[i].walker[j].pars==NULL) {
          fprintf(stderr,"Could not allocate struct pars\n");
          exit(EXIT_FAILURE);
      }
      self->ball_2[i].walker[j].pars=calloc(npars,sizeof(double));
      if (self->ball_2[i].walker[j].pars==NULL) {
          fprintf(stderr,"Could not allocate struct pars\n");
          exit(EXIT_FAILURE);
      }
    }
  }
  self->nsteps=nsteps;
  return self;
}

void free_chain(chain *c){
  int nsteps=c->nsteps;
  for(int i=0; i<nsteps; i++){
    int nwalkers=c->ball_1[i].nwalkers;
    for(int j=0; j<nwalkers; j++){
      free(c->ball_1[i].walker[j].pars);
      free(c->ball_2[i].walker[j].pars);
    }
    free(c->ball_1[i].walker);
    free(c->ball_2[i].walker);
  }
  free(c->ball_1);
  free(c->ball_2);
  free(c);
}


int main( int argc, char ** argv )
{
  int nwalkers=5;
  int npars=3;
  int nsteps=10;

  chain *my_chain=allocate_chain(nsteps,nwalkers,npars);
  free_chain(my_chain);
  return 0;
}