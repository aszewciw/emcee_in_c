#include <stdio.h>
#include <stdlib.h>
#include <stddef.h> // offsetof()
#include "mpi.h"


// this is the custom structure that will be communicated in MPI
typedef struct mca_step {
  int *accept;      // 0 or 1; [nwalkers] values
  double *lnprob;   // ln(prob) at point in param space; [nwalkers] values
  double *pars;     // values of parameters; npars*nwalkers values
} mca_step;

typedef struct mca_chain {
  int nwalkers;
  int steps_per_walker;
  int npars;
  mca_step *steps;      // one "step" is the params for the ensemble of walkers
} mca_chain;

mca_chain* allocate_chain(int nwalkers, int nsteps, int npars){
  struct mca_chain *self=calloc(1,sizeof(mca_chain));
  if (self==NULL) {
      fprintf(stderr,"Could not allocate struct mca_chain\n");
      exit(EXIT_FAILURE);
  }

  self->steps=calloc(nsteps,sizeof(mca_step));
  if (self->steps==NULL) {
      fprintf(stderr,"Could not allocate mca_chain steps\n");
      exit(EXIT_FAILURE);
  }

  for(int i=0; i<nsteps; i++){

    self->steps[i].pars=calloc(nwalkers*npars,sizeof(double));
    if (self->steps[i].pars==NULL) {
        fprintf(stderr,"Could not allocate mca_chain pars\n");
        exit(EXIT_FAILURE);
    }
    self->steps[i].lnprob=calloc(nwalkers,sizeof(double));
    if (self->steps[i].lnprob==NULL) {
        fprintf(stderr,"Could not allocate mca_chain lnprob\n");
        exit(EXIT_FAILURE);
    }
    self->steps[i].accept=calloc(nwalkers,sizeof(int));
    if (self->steps[i].accept==NULL) {
        fprintf(stderr,"Could not allocate mca_chain accept\n");
        exit(EXIT_FAILURE);
    }
  }

  self->nwalkers=nwalkers;
  self->steps_per_walker=nsteps;
  self->npars=npars;

  return self;
}

void free_chain(mca_chain *chain){
  int nsteps=chain->steps_per_walker;
  for(int i=0; i<nsteps; i++){
    free(chain->steps[i].pars);
    free(chain->steps[i].lnprob);
    free(chain->steps[i].accept);
  }
  free(chain);
}

// void printStudent( student * student_x , int prank, int np )
// {

// }

int main( int argc, char ** argv )
{
   int nprocs,rank;

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   int nwalkers=1000;
   int nsteps=500000;
   int npars=10;

   mca_chain *chain = allocate_chain(nwalkers,nsteps,npars);
   free_chain(chain);
   MPI_Finalize();
   return 0;
}