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
// ensemble* allocate_ensemble(size_t nwalkers, size_t npars){

//     struct ensemble *self=calloc(1,sizeof(ensemble));
//     if (self==NULL) {
//         fprintf(stderr,"Could not allocate struct ensemble\n");
//         exit(EXIT_FAILURE);
//     }

//     self->nwalkers=nwalkers;
//     self->npars=npars;

//     self->walker=calloc(nwalkers,sizeof(walker_pos));
//     if (self->walker==NULL) {
//         fprintf(stderr,"Could not allocate struct walker_pos\n");
//         exit(EXIT_FAILURE);
//     }
//     return self;
// }

/* -------------------------------------------------------------------------- */
void free_walkers(walker_pos *w){
    free(w);
}

/* -------------------------------------------------------------------------- */
// void free_ensemble(ensemble *e){
//     free(e->walker);
//     free(e);
// }

/* -------------------------------------------------------------------------- */
// void free_chain(chain *c){
//     size_t nsteps, i;
//     nsteps=c->nsteps;

//     for(i=0; i<nsteps; i++){
//         free(c->ensemble_A[i].walker);
//         free(c->ensemble_B[i].walker);
//     }
//     free(c->ensemble_A);
//     free(c->ensemble_B);
//     free(c);
// }

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
int walker_accept(double lnprob_old, double lnprob_new, size_t npars, double z)
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
// void step_walkers(walker_pos *walkers, ensemble *comp_walkers, size_t nwalkers,
//                   double a, double (*lnprob)(const double *, size_t, const void *),
//                   const void *userdata)
// {
//     size_t iwalker,ipar,icomp,npars,ncomp;
//     int accept;
//     double par_old, par_comp;
//     double pars_new[NPARS];
//     double lnprob_old,lnprob_new,z;

//     ncomp=comp_walkers->nwalkers;
//     npars=comp_walkers->npars;

//     for(iwalker=0; iwalker<nwalkers; iwalker++){
//         lnprob_old = walkers[iwalker].lnprob;
//         icomp = rand_walker(ncomp);
//         z = rand_gofz(a);

//         for(ipar=0; ipar<npars; ipar++){
//             par_old = walkers[iwalker].pars[ipar];
//             par_comp = comp_walkers->walker[icomp].pars[ipar];
//             pars_new[ipar] = par_comp - z*(par_comp-par_old);
//         }

//         lnprob_new = lnprob(pars_new,npars,userdata);
//         accept = walker_accept(lnprob_old,lnprob_new,npars,z);
//         walkers[iwalker].accept = accept;
//         if(accept){
//             walkers[iwalker].lnprob=lnprob_new;
//             for(ipar=0; ipar<npars; ipar++){
//                 walkers[iwalker].pars[ipar]=pars_new[ipar];
//             }
//         }
//     }
// }
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
    fprintf(file, "# step accept lnprob [pars]\n");

    fclose(file);
    return 1;
}

/* -------------------------------------------------------------------------- */
int write_step(const char *fname, const struct walker_pos *ensemble_A,
               const struct walker_pos *ensemble_B, size_t nwalkers, size_t istep)
{
    size_t npars, iwalker, ipar;

    npars=(size_t)NPARS;

    FILE *file=fopen(fname,"a");
    if (file==NULL) {
        fprintf(stderr,"Error: could not open file '%s'\n", fname);
        return 0;
    }

    for(iwalker=0; iwalker<nwalkers; iwalker++){
        fprintf(file,"%zu\t%d\t%.16g",
                istep,
                ensemble_A[iwalker].accept,
                ensemble_A[iwalker].lnprob);
        for(ipar=0; ipar<npars; ipar++){
            fprintf(file,"\t%.16g", ensemble_A[iwalker].pars[ipar]);
        }
        fprintf(file,"\n");
    }
    for(iwalker=0; iwalker<nwalkers; iwalker++){
        fprintf(file,"%zu\t%d\t%.16g",
                istep,
                ensemble_B[iwalker].accept,
                ensemble_B[iwalker].lnprob);
        for(ipar=0; ipar<npars; ipar++){
            fprintf(file,"\t%.16g", ensemble_B[iwalker].pars[ipar]);
        }
        fprintf(file,"\n");
    }

    fclose(file);
    return 1;
}

/* -------------------------------------------------------------------------- */
void create_trials(walker_pos *trial, const struct walker_pos *walkers,
                   const struct walker_pos *comp_walkers, double *z_array,
                   size_t nwalkers, double a)
{
    size_t iwalker, ipar, icomp, npars;
    double z, par_old, par_comp, par_new;

    npars=(size_t)NPARS;

    for(iwalker=0; iwalker<nwalkers; iwalker++){
        icomp = rand_walker(nwalkers);
        z = rand_gofz(a);
        z_array[iwalker] = z;

        for(ipar=0; ipar<npars; ipar++){
            par_old = walkers[iwalker].pars[ipar];
            par_comp = comp_walkers[icomp].pars[ipar];
            par_new = par_comp - z*(par_comp - par_old);
            trial[iwalker].pars[ipar] = par_new;
        }
    }
}
/* -------------------------------------------------------------------------- */
void update_positions(walker_pos *walkers, const struct walker_pos *trial,
                      double *z_array, size_t nwalkers)
{
    size_t iwalker, ipar, npars;
    double z, lnprob_old, lnprob_new;
    int accept;
    npars=(size_t)NPARS;

    for(iwalker=0; iwalker<nwalkers; iwalker++){
        z = z_array[iwalker];
        lnprob_old = walkers[iwalker].lnprob;
        lnprob_new = trial[iwalker].lnprob;
        accept = walker_accept(lnprob_old, lnprob_new, npars, z);
        walkers[iwalker].accept = accept;
        if(accept){
            walkers[iwalker].lnprob = lnprob_new;
            for(ipar=0; ipar<npars; ipar++){
                walkers[iwalker].pars[ipar] = trial[iwalker].pars[ipar];
            }
        }
    }
}
/* -------------------------------------------------------------------------- */
void manager(walker_pos *start_pos, double a, const char *fname, int nburn){

    walker_pos *ensemble_A, *ensemble_B, *trial;
    size_t nsteps, nwalkers, npars, nwalkers_over_two, iwalker, ipar, istep;
    int nprocs, rank, irecv;
    double lnprob_tmp;
    time_t t;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    int current_task[nprocs];

    nsteps   = (size_t)NSTEPS;
    nwalkers = (size_t)NWALKERS;
    npars    = (size_t)NPARS;
    nwalkers_over_two = nwalkers/2;
    double z_array[nwalkers_over_two];

    ensemble_A = allocate_walkers(nwalkers_over_two);
    ensemble_B = allocate_walkers(nwalkers_over_two);
    trial = allocate_walkers(nwalkers_over_two);

    srand((unsigned)(time(&t)));

    /* get initial lnprob values */
    iwalker=0;
    for(rank=1; rank<nprocs; rank++){
        current_task[rank]=iwalker;
        MPI_Send(&start_pos[iwalker].pars[0], npars, MPI_DOUBLE, rank, WORKTAG, MPI_COMM_WORLD);
        iwalker++;
    }

    /* receive lnprob's as they come in, and send out remaining walker positions */
    while(iwalker<nwalkers){
        MPI_Recv(&lnprob_tmp, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        irecv=current_task[status.MPI_SOURCE];
        start_pos[irecv].lnprob = lnprob_tmp;
        current_task[status.MPI_SOURCE]=iwalker;
        MPI_Send(&start_pos[iwalker].pars[0], npars, MPI_DOUBLE, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
        iwalker++;
    }

    /* receive remaining lnprobs -- should be one for each worker, but they arrive in unknown order */
    for (rank = 1; rank < nprocs; rank++) {
        MPI_Recv(&lnprob_tmp, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        irecv=current_task[status.MPI_SOURCE];
        start_pos[irecv].lnprob = lnprob_tmp;
    }

    /* fill the two ensembles with starting positions */
    for(iwalker=0; iwalker<nwalkers_over_two; iwalker++){
        ensemble_A[iwalker].accept=1;
        ensemble_B[iwalker].accept=1;
        ensemble_A[iwalker].lnprob=start_pos[iwalker].lnprob;
        ensemble_B[iwalker].lnprob=start_pos[iwalker+nwalkers_over_two].lnprob;
        for(ipar=0; ipar<npars; ipar++){
            ensemble_A[iwalker].pars[ipar]=start_pos[iwalker].pars[ipar];
            ensemble_B[iwalker].pars[ipar]=start_pos[iwalker+nwalkers_over_two].pars[ipar];
        }
    }

    write_header(fname,nsteps, nwalkers, npars);
    write_step(fname, ensemble_A, ensemble_B, nwalkers_over_two, 0);

    /* begin the chain */
    for(istep=1; istep<nsteps; istep++){
        /*----------------------------- ensemble A ---------------------------*/
        create_trials(trial, ensemble_A, ensemble_B, z_array, nwalkers_over_two, a);

        /* send out first batch of trial positions */
        iwalker=0;
        for(rank=1; rank<nprocs; rank++){
            current_task[rank]=iwalker;
            MPI_Send(&trial[iwalker].pars[0], npars, MPI_DOUBLE, rank, WORKTAG, MPI_COMM_WORLD);
            iwalker++;
        }

        /* receive lnprob's as they come in, and send out remaining walker positions */
        while(iwalker<nwalkers_over_two){
            MPI_Recv(&lnprob_tmp, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            irecv=current_task[status.MPI_SOURCE];
            trial[irecv].lnprob = lnprob_tmp;
            current_task[status.MPI_SOURCE]=iwalker;
            MPI_Send(&trial[iwalker].pars[0], npars, MPI_DOUBLE, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
            iwalker++;
        }

        /* receive remaining lnprobs -- should be one for each worker, but they arrive in unknown order */
        for (rank = 1; rank < nprocs; rank++) {
            MPI_Recv(&lnprob_tmp, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            irecv=current_task[status.MPI_SOURCE];
            trial[irecv].lnprob = lnprob_tmp;
        }

        /* update walker positions */
        update_positions(ensemble_A, trial, z_array, nwalkers_over_two);

        /*----------------------------- ensemble B ---------------------------*/
        create_trials(trial, ensemble_B, ensemble_A, z_array, nwalkers_over_two, a);

        /* send out first batch of trial positions */
        iwalker=0;
        for(rank=1; rank<nprocs; rank++){
            current_task[rank]=iwalker;
            MPI_Send(&trial[iwalker].pars[0], npars, MPI_DOUBLE, rank, WORKTAG, MPI_COMM_WORLD);
            iwalker++;
        }

        /* receive lnprob's as they come in, and send out remaining walker positions */
        while(iwalker<nwalkers_over_two){
            MPI_Recv(&lnprob_tmp, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            irecv=current_task[status.MPI_SOURCE];
            trial[irecv].lnprob = lnprob_tmp;
            current_task[status.MPI_SOURCE]=iwalker;
            MPI_Send(&trial[iwalker].pars[0], npars, MPI_DOUBLE, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
            iwalker++;
        }

        /* receive remaining lnprobs -- should be one for each worker, but they arrive in unknown order */
        for (rank = 1; rank < nprocs; rank++) {
            MPI_Recv(&lnprob_tmp, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            irecv=current_task[status.MPI_SOURCE];
            trial[irecv].lnprob = lnprob_tmp;
        }

        /* update walker positions */
        update_positions(ensemble_B, trial, z_array, nwalkers_over_two);

        write_step(fname, ensemble_A, ensemble_B, nwalkers_over_two, istep);
    }

    /* chain complete -- tell workers to exit */
    for (rank = 1; rank < nprocs; rank++) {
        MPI_Send(0, 0, MPI_DOUBLE, rank, DIETAG, MPI_COMM_WORLD);
    }

    free_walkers(ensemble_A);
    free_walkers(ensemble_B);
    free_walkers(trial);
}

/* -------------------------------------------------------------------------- */
void worker(const void *userdata, double (*lnprob)(const double *, size_t, const void *)){

    size_t npars;
    MPI_Status status;
    double pars[npars];
    double lnprob_new;

    npars=(size_t)NPARS;

    while(1){
        MPI_Recv(&pars[0], npars, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == DIETAG) {
            break;
        }
        lnprob_new = lnprob(pars, npars, userdata);
        MPI_Send(&lnprob_new, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

}
/* -------------------------------------------------------------------------- */
void run_chain_loadbalancing(int *argc, char ***argv, walker_pos *start_pos, double a,
                             double (*lnprob)(const double *, size_t, const void *),
                             const void *userdata, const char *fname, int nburn)
{
    /*========================================================================*/
    /* MPI stuff */
    int nprocs, rank, nworkers;
    MPI_Init(argc,argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* chain things */
    size_t nsteps, nwalkers, npars, nwalkers_over_two;
    double *pars;

    /*========================================================================*/
    /* check that nwalkers is even */
    nwalkers = (size_t)NWALKERS;
    if(nwalkers%2==1){
        if(rank==0){
            fprintf(stderr, "Use an even number of walkers. They will be split into two ensembles\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    nwalkers_over_two = nwalkers/2;
    nworkers = nprocs-1;

    if(nworkers>nwalkers_over_two){
        if(rank==0){
            fprintf(stderr, "Attempting to split a 'half-ensemble' of %zu walkers across %d workers.\n",
                    nwalkers_over_two, nworkers);
            fprintf(stderr, "I don't know how to do that. Change the values in the 'pars.h' file\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    if(rank==0){
        manager(start_pos, a, fname, nburn);
    }
    else{
        worker(userdata, lnprob);
    }

    /* end MPI */
    MPI_Finalize();
}

/* -------------------------------------------------------------------------- */
void run_chain(int *argc, char ***argv, walker_pos *start_pos, double a,
               double (*lnprob)(const double *, size_t, const void *),
               const void *userdata, const char *fname, int nburn,
               int load_balancing)
{
    if(load_balancing){
        run_chain_loadbalancing(argc,argv,start_pos,a,lnprob,userdata,fname,nburn);
    }
    else{
        // run_chain_standard(argc,argv,start_pos,a,lnprob,userdata,fname,nburn);
        fprintf(stderr, "Error: load_balancing must be enabled currently\n");
        exit(EXIT_FAILURE);
    }
}