// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
// #include <string.h>
// #include "randn.h"
#include "emcee.h"

/* -------------------------------------------------------------------------- */
walker_pos* allocate_walkers(int nwalkers){
    struct walker_pos *self=calloc(nwalkers,sizeof(walker_pos));
    if (self==NULL) {
        fprintf(stderr,"Could not allocate struct walker_pos\n");
        exit(EXIT_FAILURE);
    }
    return self;
}

/* -------------------------------------------------------------------------- */
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
    return self;
}

/* -------------------------------------------------------------------------- */
chain* allocate_chain(int nsteps, int nwalkers, int npars){

    int nwalkers_over_two=nwalkers/2;
    int i;
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

    for(i=0; i<nsteps; i++){
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
    }
    self->nsteps=nsteps;
    return self;
}

/* -------------------------------------------------------------------------- */
void free_walkers(walker_pos *w){
    free(w);
}

/* -------------------------------------------------------------------------- */
void free_ensemble(ensemble *e){
    free(e->walker);
    free(e);
}

/* -------------------------------------------------------------------------- */
void free_chain(chain *c){
    int nsteps=c->nsteps;
    for(int i=0; i<nsteps; i++){
        free(c->ball_1[i].walker);
        free(c->ball_2[i].walker);
    }
    free(c->ball_1);
    free(c->ball_2);
    free(c);
}

/* -------------------------------------------------------------------------- */