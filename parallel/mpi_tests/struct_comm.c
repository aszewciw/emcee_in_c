#include <stdio.h>
#include <stdlib.h>
#include <stddef.h> // offsetof()
#include "mpi.h"


// this is the custom structure that will be communicated in MPI
typedef struct {
  int *accept;      // 0 or 1; [nwalkers] values
  double *lnprob;   // ln(prob) at point in param space; [nwalkers] values
  double *pars;     // values of parameters; npars*nwalkers values
} mca_step;

typedef struct {
  int nwalkers;
  int steps_per_walker;
  int npars;
  mca_step *steps;      // one "step" is the params for the ensemble of walkers
} mca_chain;

mca_chain allocate_chain(int nwalkers, int nsteps, int npars){
  struct mca_chain *self=calloc(1,sizeof(mca_chain));
  if (self==NULL) {
      fprintf(stderr,"Could not allocate struct mca_chain\n");
      exit(EXIT_FAILURE);
  }

  self.nwalkers=nwalkers;
  self.steps_per_walker=nsteps;
  self.npars=npars;

  self->steps=calloc(nsteps,sizeof(mca_step));
  if (self->steps==NULL) {
      fprintf(stderr,"Could not allocate mca_chain steps\n");
      exit(EXIT_FAILURE);
  }

  for(int i=0; i<nsteps; i++){

    self[i]->pars=calloc(nwalkers*npars,sizeof(double));
    if (self[i]->pars==NULL) {
        fprintf(stderr,"Could not allocate mca_chain pars\n");
        exit(EXIT_FAILURE);
    }
    self[i]->lnprob=calloc(nwalkers,sizeof(double));
    if (self[i]->lnprob==NULL) {
        fprintf(stderr,"Could not allocate mca_chain lnprob\n");
        exit(EXIT_FAILURE);
    }
    self[i]->accept=calloc(nwalkers,sizeof(int));
    if (self[i]->accept==NULL) {
        fprintf(stderr,"Could not allocate mca_chain accept\n");
        exit(EXIT_FAILURE);
    }
  }

  return self;
}

void free_chain(mca_chain *chain){
  int nsteps=chain->nsteps;
  for(int i=0; i<nsteps; i++){
    free(chain[i]->pars);
    free(chain[i]->lnprob);
    free(chain[i]->accept);
  }
  free(chain);
}

void printStudent( student * student_x , int prank, int np )
{

}

int main( int argc, char ** argv )
{
   int nprocs,rank;

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   int nwalkers=1000;
   int nsteps=500000;
   int npars=10;

   if(rank==0) mca_chain *chain = allocate_chain(nwalkers,nsteps,npars);
   if(rank==0) free_chain(chain);
   MPI_Finalize();
   return 0;
}