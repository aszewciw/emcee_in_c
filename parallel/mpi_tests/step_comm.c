#include <stdio.h>
#include <stdlib.h>
#include <stddef.h> // offsetof()
#include "mpi.h"

// #define NPARS 2

// this is the custom structure that will be communicated in MPI
typedef struct walker_pos {
  int accept;
  double lnprob;
  double *pars;
  // double pars[NPARS];
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

walker_pos* allocate_walkers(int nwalkers, int npars){
  struct walker_pos *self=calloc(nwalkers,sizeof(walker_pos));
  if (self==NULL) {
    fprintf(stderr,"Could not allocate struct walker_pos\n");
    exit(EXIT_FAILURE);
  }

  for(int i=0; i<nwalkers; i++){
    self[i].pars=calloc(npars,sizeof(double));;
    if (self[i].pars==NULL) {
      fprintf(stderr,"Could not allocate array pars\n");
      exit(EXIT_FAILURE);
    }
  }
  return self;
}

ensemble* allocate_ensemble(int nwalkers, int npars){

  struct ensemble *self=calloc(1,sizeof(ensemble));
  if (self==NULL) {
    fprintf(stderr,"Could not allocate struct ensemble\n");
    exit(EXIT_FAILURE);
  }

  self->nwalkers=nwalkers;
  self->npars=npars;

  self->walker=calloc(nwalkers,sizeof(walker_pos));
  if (self->walker==NULL) {
    fprintf(stderr,"Could not allocate struct walker_pos\n");
    exit(EXIT_FAILURE);
  }

  for(int i=0; i<nwalkers; i++){
    self->walker[i].pars=calloc(npars,sizeof(double));
    if (self->walker[i].pars==NULL) {
        fprintf(stderr,"Could not allocate struct pars\n");
        exit(EXIT_FAILURE);
    }
  }
  return self;
}

chain* allocate_chain(int nsteps, int nwalkers, int npars){

  int nwalkers_over_two=nwalkers/2;
  struct chain *self=calloc(1,sizeof(chain));
  if (self==NULL) {
      fprintf(stderr,"Could not allocate struct chain\n");
      exit(EXIT_FAILURE);
  }

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
          fprintf(stderr,"Could not allocate array pars\n");
          exit(EXIT_FAILURE);
      }
      self->ball_2[i].walker[j].pars=calloc(npars,sizeof(double));
      if (self->ball_2[i].walker[j].pars==NULL) {
          fprintf(stderr,"Could not allocate array pars\n");
          exit(EXIT_FAILURE);
      }
    }
  }
  self->nsteps=nsteps;
  return self;
}

void free_walkers(walker_pos *w, int nwalkers){
  for(int i=0; i<nwalkers; i++){
    free(w[i].pars);
  }
  free(w);
}


void free_ensemble(ensemble *e){
  int nwalkers=e->nwalkers;
  for(int i=0; i<nwalkers; i++){
    free(e->walker[i].pars);
  }
  free(e->walker);
  free(e);
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
  int nprocs,rank;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int nsteps=5;
  int nwalkers=6;
  int npars=2;
  int nwalkers_over_two=nwalkers/2;

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

  chain *my_chain;
  if (rank==0){
    my_chain=allocate_chain(nsteps,nwalkers,npars);
  }

  ensemble *my_ensemble=allocate_ensemble(nwalkers,npars);
  walker_pos *my_walkers=allocate_walkers(slice_length,npars);

  // have each chain fill its walker_pos
  for(int i=0; i<slice_length; i++){
    my_walkers[i].accept = rank+1;
    my_walkers[i].lnprob = (rank+1)*10.0;
    for(int j=0; j<npars; j++){
      my_walkers[i].pars[j] = (double)(rank+i+j+1);
    }
  }

  MPI_Datatype MPI_WALKER;
  MPI_Datatype type[3] = { MPI_INT, MPI_DOUBLE, MPI_DOUBLE };
  int blocklen[3] = { 1, 1, npars }; // size of each data element in struct

  MPI_Aint disp[3]; // array of displacements; one for each data member
  disp[0] = offsetof(walker_pos,accept);
  disp[1] = offsetof(walker_pos,lnprob);
  disp[2] = offsetof(walker_pos,pars);

  MPI_Type_create_struct(3,blocklen,disp,type,&MPI_WALKER);
  MPI_Type_commit(&MPI_WALKER);

  // for(int istep=0; istep<nsteps; istep++){
  //   for(int iwalker=0; iwalker<nwalkers_over_two; iwalker++){
  //     my_chain->ball_1[istep].walker[iwalker].accept=istep*iwalker;
  //     my_chain->ball_1[istep].walker[iwalker].lnprob=(double)(istep*iwalker+100);
  //     my_chain->ball_2[istep].walker[iwalker].accept=istep*iwalker*100;
  //     my_chain->ball_2[istep].walker[iwalker].lnprob=(double)(istep*iwalker+200);
  //     for(int ipar=0;ipar<npars;ipar++){
  //       my_chain->ball_1[istep].walker[iwalker].pars[ipar]=(double)(istep*iwalker*ipar+1000);
  //       my_chain->ball_2[istep].walker[iwalker].pars[ipar]=(double)(istep*iwalker*ipar+2000);
  //     }
  //   }
  // }

  MPI_Allgatherv(&my_walkers, slice_length, MPI_WALKER,
                 &my_ensemble, slice_length, lower_ind,
                 MPI_WALKER, MPI_COMM_WORLD);

  if(rank==0) free_chain(my_chain);
  free_ensemble(my_ensemble);
  free_walkers(my_walkers,slice_length);

  MPI_Type_free(&MPI_WALKER);
  MPI_Finalize();
  return 0;
}