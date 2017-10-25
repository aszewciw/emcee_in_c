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

mca_step* allocate_step(int nwalkers, int npars){

  struct mca_step *self=calloc(1,sizeof(mca_step));
  if (self==NULL) {
      fprintf(stderr,"Could not allocate step mca_step\n");
      exit(EXIT_FAILURE);
  }

  self->pars=calloc(nwalkers*npars,sizeof(double));
  if (self->pars==NULL) {
      fprintf(stderr,"Could not allocate step pars\n");
      exit(EXIT_FAILURE);
  }
  self->lnprob=calloc(nwalkers,sizeof(double));
  if (self->lnprob==NULL) {
      fprintf(stderr,"Could not allocate step lnprob\n");
      exit(EXIT_FAILURE);
  }
  self->accept=calloc(nwalkers,sizeof(int));
  if (self->accept==NULL) {
      fprintf(stderr,"Could not allocate step accept\n");
      exit(EXIT_FAILURE);
  }

  return self;
}


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

void free_step(mca_step *step){
  free(step->pars);
  free(step->lnprob);
  free(step->accept);
  free(step);
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

  int nwalkers=6;
  int nsteps=5;
  int npars=2;

  /* Establish slice of walkers for each process to handle */
  int slice_length, lower_ind, upper_ind;
  int remain = nwalkers % nprocs;

  /* Make slices as even as possible */
  slice_length = nwalkers / nprocs;
  lower_ind = rank * slice_length;
  if (rank < remain){
      lower_ind += rank;
      slice_length++;
  }
  else lower_ind += remain;
  upper_ind = lower_ind + slice_length;

  // allocate space where we'll store full chain; rank 0 only
  // idk if I even really want to store the whole chain. Might be good for quickly doing some stats.
  /*
  It might actually be good to leave one proc free for i/o and other meta tasks.
  For now, I'm not going to do this.
  */
  mca_chain *chain;
  if (rank==0){
    chain = allocate_chain(nwalkers,nsteps,npars);
  }

  /*
  This will store one step for nwalkers; Each proc has its own copy, but is only
  responsible for its slice of walkers; i.e., the rest of this structure will need
  to be communicated via MPI
  */
  mca_step *step = allocate_step(nwalkers,npars);
  mca_step *step_combined = allocate_step(nwalkers,npars);

  // have each chain fill its one step
  for(int i=lower_ind; i<upper_ind; i++){
    step->accept[i] = rank;
    step->lnprob[i] = rank*10.0;
    for(int j=0; j<npars; j++){
      step->pars[i+j] = (double)(rank+i+j);
    }
  }

  MPI_Datatype MCA_STEP;
  MPI_Datatype type[3] = { MPI_INT, MPI_DOUBLE, MPI_DOUBLE };
  int blocklen[3] = { nwalkers, nwalkers, nwalkers*npars }; // size of each data element in struct

  MPI_Aint disp[3]; // array of displacements; one for each data member
  disp[0] = offsetof(mca_step,accept);
  disp[1] = offsetof(mca_step,lnprob);
  disp[2] = offsetof(mca_step,pars);

  MPI_Type_create_struct(3,blocklen,disp,type,&MCA_STEP);
  MPI_Type_commit(&MCA_STEP);

  MPI_Allreduce(step,step_combined,1,MCA_STEP,MPI_SUM,MPI_COMM_WORLD);

  int current_rank=0;
  while(current_rank<nprocs){
    if (current_rank==rank){
      fprintf(stderr, "Rank %d is printing:\n\n", current_rank);
      for(int i=0; i<nwalkers; i++){
        fprintf(stderr, "Walker %d accept: %d\n", i,step->accept[i]);
        fprintf(stderr, "Walker %d lnprob: %lf\n", i,step->lnprob[i]);
        for(int j=0; j<npars; j++){
          fprintf(stderr, "Walker %d param %d: %lf\n", i,j,step->pars[i+j]);
        }
      }
    }
    current_rank+=1;
    MPI_Barrier(MPI_COMM_WORLD);
  }

  if (rank==0) {
    free_chain(chain);
  }
  free_step(step);
  free_step(step_combined);

  MPI_Type_free(&MCA_STEP);
  MPI_Finalize();
  return 0;
}