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
    double r = rand_0to1();

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

/* -------------------------------------------------------------------------- */
int write_chain(const struct chain *c, const char *fname)
{
    int nwalkers, nsteps, npars, nwalkers_over_two;
    int istep,iwalker,ipar;

    nsteps = c->nsteps;
    nwalkers_over_two = c->ball_1[0].nwalkers;
    nwalkers = nwalkers_over_two*2;
    npars = c->ball_1[0].npars;


    FILE *file=fopen(fname,"w");
    if (file==NULL) {
        fprintf(stderr,"Error: could not open file '%s'\n", fname);
        return 0;
    }

    fprintf(file, "# nsteps\tnwalkers\tnpars\n");
    fprintf(file,"%d\t%d\t%d\n", nsteps, nwalkers, npars);
    fprintf(file, "# accept\tlnprob\tpars\n");
    for (istep=0; istep<nsteps; istep++) {
        for(iwalker=0; iwalker<nwalkers_over_two; iwalker++){
            fprintf(file,"%d\t%.16g\t",
                    c->ball_1[istep].walker[iwalker].accept,
                    c->ball_1[istep].walker[iwalker].lnprob);
            for(ipar=0; ipar<npars; ipar++){
                fprintf(file,"%.16g\t",
                        c->ball_1[istep].walker[iwalker].pars[ipar]);
            }
            fprintf(file,"\n");
        }
        for(iwalker=0; iwalker<nwalkers_over_two; iwalker++){
            fprintf(file,"%d\t%.16g\t",
                    c->ball_2[istep].walker[iwalker].accept,
                    c->ball_2[istep].walker[iwalker].lnprob);
            for(ipar=0; ipar<npars; ipar++){
                fprintf(file,"%.16g\t",
                        c->ball_2[istep].walker[iwalker].pars[ipar]);
            }
            fprintf(file,"\n");
        }
    }

    fclose(file);
    return 1;
}


void run_chain(int *argc, char ***argv, double a, double *centers, double *widths,
               double (*lnprob)(const double *, int, const void *),
               const void *userdata, const char *fname)
{
    /*========================================================================*/
    /* MPI stuff */
    int nprocs, rank;
    MPI_Init(argc,argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* generic indices */
    int iproc, iwalker, ipar, istep, iensemble;

    /* chain things */
    int nsteps, nwalkers, npars, nwalkers_over_two;
    double *pars;
    ensemble *ensemble_A, *ensemble_B;
    walker_pos *my_walkers, *start_pos;
    chain *my_chain;

    /* more MPI stuff */
    int slice_length, lower_ind, upper_ind, remain;
    int mpi_disp[nprocs], counts[nprocs];
    int tmp_start, tmp_slice;

    /* rng */
    time_t t;

    /*========================================================================*/
    /* Set up basic chain parameters */
    nsteps   = (int)NSTEPS;
    nwalkers = (int)NWALKERS;
    npars    = (int)NPARS;
    if(nwalkers%2==1) nwalkers++;   // make nwalkers even

    nwalkers_over_two = nwalkers/2;

    if(nprocs>nwalkers_over_two){
        if(rank==0){
            fprintf(stderr, "Attempting to split a 'half-ensemble' of %d walkers across %d processes is silly.\n",
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

    /* set up chain, ensembles, walkers */
    if (rank==0) my_chain = allocate_chain(nsteps,nwalkers,npars);

    ensemble_A = allocate_ensemble(nwalkers_over_two,npars);
    ensemble_B = allocate_ensemble(nwalkers_over_two,npars);

    my_walkers = allocate_walkers(slice_length);

    /*========================================================================*/

    /* Get things ready to begin the chain */
    /* set up the rng here */
    srand((unsigned)(time(&t))+rank);

    /* Have each process make its own guess. We'll overwrite it in a second */
    start_pos = make_guess(centers,widths,nwalkers,npars);

    /* Have process 0 send data to all others */
    MPI_Bcast(&start_pos[0], nwalkers, MPI_WALKER, 0, MPI_COMM_WORLD);

    /* Fill the two ensembles with walker positions */
    for(iwalker=0; iwalker<nwalkers_over_two; iwalker++){
        for(ipar=0; ipar<npars; ipar++){
            ensemble_A->walker[iwalker].pars[ipar]=start_pos[iwalker].pars[ipar];
            ensemble_B->walker[iwalker].pars[ipar]=start_pos[iwalker+nwalkers_over_two].pars[ipar];
        }
    }

    /* Fill each proc's walkers with appropriate starting positions */
    /* We'll first do ensemble_B */
    /* also calculate lnprob */
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

    /* fill first step of chain */
    if(rank==0){
        for(iwalker=0;iwalker<nwalkers_over_two;iwalker++){
            my_chain->ball_1[0].walker[iwalker].accept=ensemble_A->walker[iwalker].accept;
            my_chain->ball_1[0].walker[iwalker].lnprob=ensemble_A->walker[iwalker].lnprob;
            my_chain->ball_2[0].walker[iwalker].accept=ensemble_B->walker[iwalker].accept;
            my_chain->ball_2[0].walker[iwalker].lnprob=ensemble_B->walker[iwalker].lnprob;
            for(ipar=0;ipar<npars;ipar++){
                my_chain->ball_1[0].walker[iwalker].pars[ipar]=ensemble_A->walker[iwalker].pars[ipar];
                my_chain->ball_2[0].walker[iwalker].pars[ipar]=ensemble_B->walker[iwalker].pars[ipar];
            }
        }
    }

    /*====================== burn in would happen here =======================*/


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

        /* put A and B in chain */
        /* fill step of chain */
        if(rank==0){
            for(iwalker=0;iwalker<nwalkers_over_two;iwalker++){
                my_chain->ball_1[istep].walker[iwalker].accept=ensemble_A->walker[iwalker].accept;
                my_chain->ball_1[istep].walker[iwalker].lnprob=ensemble_A->walker[iwalker].lnprob;
                my_chain->ball_2[istep].walker[iwalker].accept=ensemble_B->walker[iwalker].accept;
                my_chain->ball_2[istep].walker[iwalker].lnprob=ensemble_B->walker[iwalker].lnprob;
                for(ipar=0;ipar<npars;ipar++){
                    my_chain->ball_1[istep].walker[iwalker].pars[ipar]=ensemble_A->walker[iwalker].pars[ipar];
                    my_chain->ball_2[istep].walker[iwalker].pars[ipar]=ensemble_B->walker[iwalker].pars[ipar];
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    /*========================================================================*/
    /* write file */
    if(rank==0) write_chain(my_chain, fname);

    MPI_Barrier(MPI_COMM_WORLD);
    /*========================================================================*/

    /* free chain */
    if(rank==0) free_chain(my_chain);
    free_ensemble(ensemble_A);
    free_ensemble(ensemble_B);
    free_walkers(my_walkers);

    /* end MPI */
    MPI_Type_free(&MPI_WALKER);
    MPI_Finalize();
}