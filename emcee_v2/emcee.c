#include "emcee.h"

/* -------------------------------------------------------------------------- */
walker_pos* allocate_walkers(size_t nwalkers){
    struct walker_pos *self=calloc(nwalkers,sizeof(walker_pos));
    if (self==NULL) {
        fprintf(stderr,"Could not allocate struct walker_pos\n");
        exit(EXIT_FAILURE);
    }
    return self;
}

/* -------------------------------------------------------------------------- */
ensemble* allocate_ensemble(size_t nwalkers, size_t npars){

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
    size_t nsteps, i;
    nsteps=c->nsteps;

    for(i=0; i<nsteps; i++){
        free(c->ensemble_A[i].walker);
        free(c->ensemble_B[i].walker);
    }
    free(c->ensemble_A);
    free(c->ensemble_B);
    free(c);
}

/* -------------------------------------------------------------------------- */
walker_pos *make_guess(double *centers, double *widths, size_t nwalkers, size_t npars)
{
    size_t ipar, iwalker;
    double center, width, val;
    walker_pos * guess=allocate_walkers(nwalkers);

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
int walker_accept(double lnprob_old,double lnprob_new,size_t npars,double z)
{
    double lnprob_diff = (npars - 1.)*log(z) + lnprob_new - lnprob_old;
    double r = rand_0to1();

    if (lnprob_diff > log(r)) {
        return 1;
    } else {
        return 0;
    }
}
/* -------------------------------------------------------------------------- */
void step_walkers(walker_pos *walkers, ensemble *comp_walkers, size_t nwalkers,
                  double a, double (*lnprob)(const double *, size_t, const void *),
                  const void *userdata)
{
    size_t iwalker,ipar,icomp,npars,ncomp;
    int accept;
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
double normal_rand()
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
size_t rand_walker(size_t n)
{
    /* fix this later for sake of uniformity */
    size_t i=rand()%n;
    return i;
}

/* -------------------------------------------------------------------------- */
double rand_gofz(double a)
{
    /*
    g(z) ~ 1/sqrt(z) on z=[1/a, a]
    After normalizing, drawing z from this distribution is done as
        z = ( (a-1) rand + 1 )^2 / a
    */

    double tmp,z;

    tmp = (a - 1.)*rand_0to1() + 1.;
    z = tmp*tmp/a;

    return z;
}

/* -------------------------------------------------------------------------- */
int write_header(const char *fname, size_t nsteps, size_t nwalkers, size_t npars)
{
    FILE *file=fopen(fname,"w");
    if (file==NULL) {
        fprintf(stderr,"Error: could not open file '%s'\n", fname);
        return 0;
    }
    fprintf(file, "# nsteps=%zu nwalkers=%zu npars=%zu\n",nsteps, nwalkers, npars);
    fprintf(file, "# accept lnprob [pars]\n");

    fclose(file);
    return 1;
}

/* -------------------------------------------------------------------------- */
int write_step(const char *fname, const struct ensemble *ensemble_A,
               const struct ensemble *ensemble_B)
{
    size_t nwalkers,npars,nwalkers_over_two;
    size_t iwalker,ipar;

    nwalkers_over_two = ensemble_A->nwalkers;
    nwalkers = nwalkers_over_two*2;
    npars = ensemble_A->npars;

    FILE *file=fopen(fname,"a");
    if (file==NULL) {
        fprintf(stderr,"Error: could not open file '%s'\n", fname);
        return 0;
    }

    for(iwalker=0; iwalker<nwalkers_over_two; iwalker++){
        fprintf(file,"%d\t%.16g\t",
                ensemble_A->walker[iwalker].accept,
                ensemble_A->walker[iwalker].lnprob);
        for(ipar=0; ipar<npars; ipar++){
            fprintf(file,"%.16g\t",
                    ensemble_A->walker[iwalker].pars[ipar]);
        }
        fprintf(file,"\n");
    }
    for(iwalker=0; iwalker<nwalkers_over_two; iwalker++){
        fprintf(file,"%d\t%.16g\t",
                ensemble_B->walker[iwalker].accept,
                ensemble_B->walker[iwalker].lnprob);
        for(ipar=0; ipar<npars; ipar++){
            fprintf(file,"%.16g\t",
                    ensemble_B->walker[iwalker].pars[ipar]);
        }
        fprintf(file,"\n");
    }

    fclose(file);
    return 1;
}

/* -------------------------------------------------------------------------- */
void run_chain(int *argc, char ***argv, walker_pos *start_pos, double a,
               double (*lnprob)(const double *, size_t, const void *),
               const void *userdata, const char *fname, int nburn)
{
    /*========================================================================*/
    /* MPI stuff */
    int nprocs, rank;
    MPI_Init(argc,argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* generic indices */
    int iproc;
    size_t iwalker, ipar, istep, iensemble;

    /* chain things */
    size_t nsteps, nwalkers, npars, nwalkers_over_two;
    double *pars;
    ensemble *ensemble_A, *ensemble_B;
    walker_pos *my_walkers;

    /* more MPI stuff */
    size_t slice_length, lower_ind, upper_ind, remain;
    int mpi_disp[nprocs], counts[nprocs];

    /* rng */
    time_t t;

    /*========================================================================*/
    /* Set up basic chain parameters */
    nsteps   = (size_t)NSTEPS;
    nwalkers = (size_t)NWALKERS;
    npars    = (size_t)NPARS;
    if(nwalkers%2==1) nwalkers++;   // make nwalkers even

    nwalkers_over_two = nwalkers/2;

    if(nprocs>nwalkers_over_two){
        if(rank==0){
            fprintf(stderr, "Attempting to split a 'half-ensemble' of %zu walkers across %d processes.\n",
                    nwalkers_over_two, nprocs);
            fprintf(stderr, "I don't know how to do that. Change the values in the 'pars.h' file\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    /*========================================================================*/
    /* make an MPI custom structure */
    int blocklen[3] = {1, 1, npars};
    MPI_Datatype MPI_WALKER;
    MPI_Datatype type[3] = { MPI_INT, MPI_DOUBLE, MPI_DOUBLE };
    MPI_Aint disp[3];

    disp[0] = offsetof(walker_pos,accept);
    disp[1] = offsetof(walker_pos,lnprob);
    disp[2] = offsetof(walker_pos,pars);
    MPI_Type_create_struct(3,blocklen,disp,type,&MPI_WALKER);
    MPI_Type_commit(&MPI_WALKER);

    /*========================================================================*/

    /* Establish slice of walkers for each process to handle */
    remain=nwalkers_over_two%nprocs;

    /* Fill mpi_disp and counts for MPI comm */
    for(iproc=0; iproc<nprocs; iproc++)
    {
        slice_length = nwalkers_over_two/nprocs;
        lower_ind = iproc*slice_length;
        if(iproc<remain){
            lower_ind+=iproc;
            slice_length++;
        }
        else{
            lower_ind+=remain;
        }
        mpi_disp[iproc] = lower_ind;
        counts[iproc] = slice_length;
    }

    /* Have each proc set its slice indices */
    lower_ind = mpi_disp[rank];
    slice_length = counts[rank];
    upper_ind = lower_ind + slice_length;

    /*========================================================================*/
    /* set up ensembles, walkers */
    ensemble_A = allocate_ensemble(nwalkers_over_two,npars);
    ensemble_B = allocate_ensemble(nwalkers_over_two,npars);

    my_walkers = allocate_walkers(slice_length);

    /*========================================================================*/

    /* Get things ready to begin the chain */
    /* set up the rng here -- ensure each proc has a different seed */
    srand((unsigned)(time(&t))+rank);

    /* Fill the two ensembles with walker positions */
    /* No need for broadcasting here because all procs have the same copy of start_pos from main */
    for(iwalker=0; iwalker<nwalkers_over_two; iwalker++){
        for(ipar=0; ipar<npars; ipar++){
            ensemble_A->walker[iwalker].pars[ipar]=start_pos[iwalker].pars[ipar];
            ensemble_B->walker[iwalker].pars[ipar]=start_pos[iwalker+nwalkers_over_two].pars[ipar];
        }
    }

    /* Fill each proc's walkers with appropriate starting positions and get lnprob */
    /* We'll first do ensemble_B */
    for(iwalker=0; iwalker<slice_length; iwalker++){
        iensemble = iwalker+lower_ind;
        for(ipar=0; ipar<npars; ipar++){
            my_walkers[iwalker].pars[ipar] = ensemble_B->walker[iensemble].pars[ipar];
        }
        my_walkers[iwalker].accept=1;
        pars=my_walkers[iwalker].pars;
        my_walkers[iwalker].lnprob=lnprob(pars,npars,userdata);
    }

    /* gather data from each proc */
    MPI_Allgatherv(&my_walkers[0], slice_length, MPI_WALKER,
                   &ensemble_B->walker[0], counts, mpi_disp,
                   MPI_WALKER, MPI_COMM_WORLD);


    /* Now we'll do ensemble_A */
    for(iwalker=0; iwalker<slice_length; iwalker++){
        iensemble = iwalker+lower_ind;
        for(ipar=0; ipar<npars; ipar++){
            my_walkers[iwalker].pars[ipar] = ensemble_A->walker[iensemble].pars[ipar];
        }
        my_walkers[iwalker].accept=1;
        pars=my_walkers[iwalker].pars;
        my_walkers[iwalker].lnprob=lnprob(pars,npars,userdata);
    }

    /* gather data from each proc */
    MPI_Allgatherv(&my_walkers[0], slice_length, MPI_WALKER,
                   &ensemble_A->walker[0], counts, mpi_disp,
                   MPI_WALKER, MPI_COMM_WORLD);

    /*=============================== burn-in ================================*/
    for(istep=0; istep<nburn; istep++){

        /* set walkers to positions for ensemble_A */
        for(iwalker=0; iwalker<slice_length; iwalker++){
            iensemble = iwalker+lower_ind;
            my_walkers[iwalker].accept = ensemble_A->walker[iensemble].accept;
            my_walkers[iwalker].lnprob = ensemble_A->walker[iensemble].lnprob;
            for(ipar=0; ipar<npars; ipar++){
                my_walkers[iwalker].pars[ipar] = ensemble_A->walker[iensemble].pars[ipar];
            }
        }
        /* move my_walkers based on ensemble_B */

        step_walkers(my_walkers, ensemble_B, slice_length, a, lnprob, userdata);

        /* allgatherv new positions of A */
        MPI_Allgatherv(&my_walkers[0], slice_length, MPI_WALKER,
                       &ensemble_A->walker[0], counts, mpi_disp,
                       MPI_WALKER, MPI_COMM_WORLD);

        /* set walkers to positions for ensemble_B */
        for(iwalker=0; iwalker<slice_length; iwalker++){
            iensemble = iwalker+lower_ind;
            my_walkers[iwalker].accept = ensemble_B->walker[iensemble].accept;
            my_walkers[iwalker].lnprob = ensemble_B->walker[iensemble].lnprob;
            for(ipar=0; ipar<npars; ipar++){
                my_walkers[iwalker].pars[ipar] = ensemble_B->walker[iensemble].pars[ipar];
            }
        }

        /* move ensemble_B based on ensemble_A */
        step_walkers(my_walkers, ensemble_A, slice_length, a, lnprob, userdata);

        /* allgatherv new positions of B */
        MPI_Allgatherv(&my_walkers[0], slice_length, MPI_WALKER,
                       &ensemble_B->walker[0], counts, mpi_disp,
                       MPI_WALKER, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
    }

    // write header and first step
    if(rank==0) write_header(fname,nsteps,nwalkers,npars);
    if(rank==0) write_step(fname,ensemble_A,ensemble_B);

    /*========================================================================*/
    /* begin the chain */
    for(istep=1; istep<nsteps; istep++){

        /* set walkers to positions for ensemble_A */
        for(iwalker=0; iwalker<slice_length; iwalker++){
            iensemble = iwalker+lower_ind;
            my_walkers[iwalker].accept = ensemble_A->walker[iensemble].accept;
            my_walkers[iwalker].lnprob = ensemble_A->walker[iensemble].lnprob;
            for(ipar=0; ipar<npars; ipar++){
                my_walkers[iwalker].pars[ipar] = ensemble_A->walker[iensemble].pars[ipar];
            }
        }
        /* move my_walkers based on ensemble_B */

        step_walkers(my_walkers, ensemble_B, slice_length, a, lnprob, userdata);

        /* allgatherv new positions of A */
        MPI_Allgatherv(&my_walkers[0], slice_length, MPI_WALKER,
                       &ensemble_A->walker[0], counts, mpi_disp,
                       MPI_WALKER, MPI_COMM_WORLD);

        /* set walkers to positions for ensemble_B */
        for(iwalker=0; iwalker<slice_length; iwalker++){
            iensemble = iwalker+lower_ind;
            my_walkers[iwalker].accept = ensemble_B->walker[iensemble].accept;
            my_walkers[iwalker].lnprob = ensemble_B->walker[iensemble].lnprob;
            for(ipar=0; ipar<npars; ipar++){
                my_walkers[iwalker].pars[ipar] = ensemble_B->walker[iensemble].pars[ipar];
            }
        }

        /* move ensemble_B based on ensemble_A */
        step_walkers(my_walkers, ensemble_A, slice_length, a, lnprob, userdata);

        /* allgatherv new positions of B */
        MPI_Allgatherv(&my_walkers[0], slice_length, MPI_WALKER,
                       &ensemble_B->walker[0], counts, mpi_disp,
                       MPI_WALKER, MPI_COMM_WORLD);

        if(rank==0) write_step(fname,ensemble_A,ensemble_B);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    free_ensemble(ensemble_A);
    free_ensemble(ensemble_B);
    free_walkers(my_walkers);

    /* end MPI */
    MPI_Type_free(&MPI_WALKER);
    MPI_Finalize();
}