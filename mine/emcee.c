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
    int i;
    for(i=0; i<nsteps; i++){
        free(c->ball_1[i].walker);
        free(c->ball_2[i].walker);
    }
    free(c->ball_1);
    free(c->ball_2);
    free(c);
}

/* -------------------------------------------------------------------------- */
struct walker_pos *make_guess(double *centers, double *widths, int nwalkers, int npars)
{
    int ipar, iwalker;
    double center, width, val;
    struct walker_pos *guess=allocate_walkers(nwalkers);

    for (ipar=0; ipar<npars; ipar++) {

        center = centers[ipar];
        width = widths[ipar];

        for (iwalker=0; iwalker<nwalkers; iwalker++) {
            val = center+width*(rand_0to1()*2.0-1.0);
            guess[iwalker].pars[ipar] = val;
        }
    }
    return guess;
}
/* -------------------------------------------------------------------------- */
// void stretch_move(double z, int npars, const double *pars_old, const double *pars_comp,
//                   double *pars_new)
// {
//     int i;
//     double val, cval;

//     for (i=0; i<npars; i++) {

//         double val=pars_old[i];
//         double cval=pars_comp[i];

//         pars_new[i] = cval - z*(cval-val);
//     }
// }
/* -------------------------------------------------------------------------- */
int walker_accept(double lnprob_old,double lnprob_new,int npars,double z)
{
    double lnprob_diff = (npars - 1.)*log(z) + lnprob_new - lnprob_old;
    double r = randu();

    if (lnprob_diff > log(r)) {
        return 1;
    } else {
        return 0;
    }
}
/* -------------------------------------------------------------------------- */
void step_walkers(walker_pos *walkers, ensemble *comp_walkers, int nwalkers,
                  double a, double (*lnprob)(const double *, int, const void *),
                  const void *userdata)
{
    int iwalker,ipar,icomp;
    int npars,ncomp,accept;
    double par_old, par_comp;
    double pars_new[NPARS];
    double lnprob_old,lnprob_new,z;

    ncomp=comp_walkers->nwalkers;
    npars=comp_walkers->npars;

    for(iwalker=0; iwalker<nwalkers; iwalker++){
        lnprob_old = walkers[iwalker].lnprob;
        icomp = rand_walker(ncomp);
        z = rand_gofz(a);

        for(ipar=0; ipar<npars; ipar++){
            par_old = walkers[iwalker].pars[ipar];
            par_comp = comp_walkers->walker[icomp].pars[ipar];
            pars_new[ipar] = par_comp - z*(par_comp-par_old);
        }

        lnprob_new = lnprob(pars_new,npars,userdata);
        accept = walker_accept(lnprob_old,lnprob_new,npars,z);
        walkers[iwalker].accept = accept;
        if(accept){
            walkers[iwalker].lnprob=lnprob_new;
            for(ipar=0; ipar<npars; ipar++){
                walkers[iwalker].pars[ipar]=pars_new[ipar];
            }
        }
    }
}
/* -------------------------------------------------------------------------- */
double rand_0to1()
{
    return (double)rand() / (double)RAND_MAX ;
}

/* -------------------------------------------------------------------------- */
double mca_randn()
{
    double x1, x2, w, y1;//, y2;

    do {
        x1 = 2.*rand_0to1() - 1.0;
        x2 = 2.*rand_0to1() - 1.0;
        w = x1*x1 + x2*x2;
    } while ( w >= 1.0 );

    w = sqrt( (-2.*log( w ) ) / w );
    y1 = x1*w;
    //y2 = x2*w;
    return y1;
}

/* -------------------------------------------------------------------------- */
int rand_walker(int n)
{
    /* fix this later for sake of uniformity */
    int i=rand()%n;
    return i;
}

/* -------------------------------------------------------------------------- */
double rand_gofz(double a)
{
    // ( (a-1) rand + 1 )^2 / a;

    double z = (a - 1.)*rand_0to1() + 1.;

    z = z*z/a;

    return z;
}