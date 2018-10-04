#include "emcee.h"

/* -------------------------------------------------------------------------- */
size_t getNlines(const char *fname, const char comment)
{
    FILE *fp= NULL;
    const int MAXLINESIZE = 10000;
    size_t nlines=0;
    char str_line[MAXLINESIZE];

    fp = fopen(fname,"rt");

    while(1){
        if(fgets(str_line, MAXLINESIZE,fp)!=NULL) {
            //WARNING: this does not remove white-space. You might
            //want to implement that (was never an issue for me)
            if(str_line[0] !=comment) nlines++;
        }
        else{
            break;
        }
    }
    fclose(fp);
    return nlines;
}

/* -------------------------------------------------------------------------- */
walker_pos* allocate_walkers(int nwalkers, int npars, int nextra)
{
    int iwalker;
    iwalker=nextra; //Just doing this to suppress warning about nextra being unused when ndef WRITE_EXTRA_DOUBLES
    struct walker_pos *self=calloc(nwalkers,sizeof(walker_pos));
    if (self==NULL) {
        fprintf(stderr,"Could not allocate struct walker_pos\n");
        exit(EXIT_FAILURE);
    }
    for(iwalker=0; iwalker<nwalkers; iwalker++){
        self[iwalker].pars = calloc(npars,sizeof(double));
        if (self[iwalker].pars==NULL) {
            fprintf(stderr,"Could not allocate pars within walker_pos\n");
            exit(EXIT_FAILURE);
        }
#ifdef WRITE_EXTRA_DOUBLES
        self[iwalker].extra_doubles = calloc(nextra,sizeof(double));
        if (self[iwalker].extra_doubles==NULL) {
            fprintf(stderr,"Could not allocate extra_doubles within walker_pos\n");
            exit(EXIT_FAILURE);
        }
#endif
    }
    return self;
}

/* -------------------------------------------------------------------------- */
void free_walkers(int nwalkers, walker_pos *w)
{
    int iwalker;
    for(iwalker=0;iwalker<nwalkers;iwalker++){
        free(w[iwalker].pars);
#ifdef WRITE_EXTRA_DOUBLES
        free(w[iwalker].extra_doubles);
#endif
    }
    free(w);
}

/* -------------------------------------------------------------------------- */
walker_pos *make_guess(int nwalkers, int npars, int nextra, double *centers,
    double *widths)
{
    int ipar, iwalker;
    double center, width, val;
    ipar=nextra;

    walker_pos * guess=allocate_walkers(nwalkers,npars,nextra);

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
int walker_accept(int npars, double lnprob_old, double lnprob_new, double z)
{
    double lnprob_diff, r;

    /* as I will recommend in documentation, return NAN for values outside priors */
    if(isnan(lnprob_new)){
        return 0;
    }

    lnprob_diff = (npars - 1.)*log(z) + lnprob_new - lnprob_old;
    r = rand_0to1();

    if (lnprob_diff > log(r)) {
        return 1;
    }
    else{
        return 0;
    }
}

/* -------------------------------------------------------------------------- */
double rand_0to1(void)
{
    // I think this is sufficient...honestly, who actually would use this with RAND_MAX~32000?
    return (double)rand() / (double)RAND_MAX ;
}

/* -------------------------------------------------------------------------- */
double normal_rand(void)
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
// taken from https://stackoverflow.com/questions/2509679/how-to-generate-a-random-integer-number-from-within-a-range
unsigned long rand_walker(unsigned long max)
{
    unsigned long num_bins, num_rand, bin_size, defect;
    long x;

    num_bins = (unsigned long) max;
    num_rand = (unsigned long) RAND_MAX + 1,
    bin_size = num_rand / num_bins,
    defect   = num_rand % num_bins;
    do
    {
        x = rand();
    }
    while (num_rand - defect <= (unsigned long)x);
    return x/bin_size;
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
int write_header(int nsteps, int nwalkers, int npars, int nextra, const char *fname)
{
    if 0{
        fprintf(stderr, "This is to suppress a warning....%d\n", nextra);
    }
    FILE *file=fopen(fname,"w");
    if (file==NULL) {
        fprintf(stderr,"Error: could not open file '%s'\n", fname);
        return 0;
    }
#ifdef WRITE_EXTRA_DOUBLES
    fprintf(file, "# nsteps=%d nwalkers=%d npars=%d\n",nsteps, nwalkers, npars);
    fprintf(file, "# step proc_id accept lnprob [pars]\n");
#else
    fprintf(file, "# nsteps=%d nwalkers=%d npars=%d nextra=%d\n",nsteps, nwalkers, npars, nextra);
    fprintf(file, "# step proc_id accept lnprob [pars] [extra_doubles]\n");
#endif
    fclose(file);
    return 1;
}

/* -------------------------------------------------------------------------- */
int write_step(int npars, int nwalkers, int istep, int nextra,
               const struct walker_pos *ensemble_A,
               const struct walker_pos *ensemble_B, const char *fname)
{
    int iwalker, ipar;
    iwalker=nextra; //Just doing this to suppress warning about nextra being unused when ndef WRITE_EXTRA_DOUBLES

    FILE *file=fopen(fname,"a");
    if (file==NULL) {
        fprintf(stderr,"Error: could not open file '%s'\n", fname);
        return 0;
    }

    for(iwalker=0; iwalker<nwalkers; iwalker++){
        fprintf(file,"%d\t%d\t%d\t%.16g",
                istep,
                ensemble_A[iwalker].rank,
                ensemble_A[iwalker].accept,
                ensemble_A[iwalker].lnprob);
        for(ipar=0; ipar<npars; ipar++){
            fprintf(file,"\t%.16g", ensemble_A[iwalker].pars[ipar]);
        }
#ifdef WRITE_EXTRA_DOUBLES
        for(int iextra=0; iextra<nextra; iextra++){
            fprintf(file, "\t%.16g", ensemble_A[iwalker].extra_doubles[iextra]);
        }
#endif
        fprintf(file,"\n");
    }
    for(iwalker=0; iwalker<nwalkers; iwalker++){
        fprintf(file,"%d\t%d\t%d\t%.16g",
                istep,
                ensemble_B[iwalker].rank,
                ensemble_B[iwalker].accept,
                ensemble_B[iwalker].lnprob);
        for(ipar=0; ipar<npars; ipar++){
            fprintf(file,"\t%.16g", ensemble_B[iwalker].pars[ipar]);
        }
#ifdef WRITE_EXTRA_DOUBLES
        for(int iextra=0; iextra<nextra; iextra++){
            fprintf(file, "\t%.16g", ensemble_B[iwalker].extra_doubles[iextra]);
        }
#endif
        fprintf(file,"\n");
    }

    fclose(file);
    return 1;
}

/* -------------------------------------------------------------------------- */
void create_trials(int nwalkers, int npars, double a, double *z_array,
                   walker_pos *trial, const struct walker_pos *walkers,
                   const struct walker_pos *comp_walkers)
{
    int iwalker, ipar, icomp;
    double z, par_old, par_comp, par_new;

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
void update_positions(int nwalkers, int npars, int nextra, double *z_array,
                      walker_pos *walkers, const struct walker_pos *trial)
{
    int iwalker, ipar, accept;
    double z, lnprob_old, lnprob_new;

    iwalker=nextra; //Just doing this to suppress warning about nextra being unused when ndef WRITE_EXTRA_DOUBLES

    for(iwalker=0; iwalker<nwalkers; iwalker++){
        z = z_array[iwalker];
        lnprob_old = walkers[iwalker].lnprob;
        lnprob_new = trial[iwalker].lnprob;
        accept = walker_accept(npars, lnprob_old, lnprob_new, z);
        walkers[iwalker].accept = accept;
        if(accept){
            walkers[iwalker].lnprob = lnprob_new;
            for(ipar=0; ipar<npars; ipar++){
                walkers[iwalker].pars[ipar] = trial[iwalker].pars[ipar];
            }
#ifdef WRITE_EXTRA_DOUBLES
            for(int iextra=0; iextra<nextra; iextra++){
                walkers[iwalker].extra_doubles[iextra]=trial[iwalker].extra_doubles[iextra];
            }
#endif
        }
    }
}

/* -------------------------------------------------------------------------- */
void manager(int nwalkers, int nsteps, int npars, int nextra, int nburn, int resume,
             double a, walker_pos *start_pos, const char *fname)
{
    int nwalkers_over_two, iwalker, ipar, istep, istart, nprocs, rank, irecv;
    const int MAXLINESIZE = 10000;
    size_t iline, Nlines, startline, offset;
    double lnprob_tmp;
    int *current_task;
    double *z_array;
    time_t t;
    char buffer[MAXLINESIZE];
    MPI_Status status;
    walker_pos *ensemble_A, *ensemble_B, *trial;
    FILE *file;

#ifdef WRITE_EXTRA_DOUBLES
    int iextra;
    double *extra_doubles = calloc(nextra+1,sizeof(double));
#endif

    iwalker=nextra; //Just doing this to suppress warning about nextra being unused when ndef WRITE_EXTRA_DOUBLES

    /* current_task let's the manager know which proc is responsible for which walker */
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    current_task = calloc(nprocs, sizeof(int));

    nwalkers_over_two = nwalkers/2;
    z_array = calloc(nwalkers_over_two, sizeof(double));

    ensemble_A = allocate_walkers(nwalkers_over_two, npars, nextra);
    ensemble_B = allocate_walkers(nwalkers_over_two, npars, nextra);
    trial = allocate_walkers(nwalkers_over_two, npars, nextra);

    srand((unsigned)(time(&t)));

    if(resume!=1){
        /* begin chain from passed start_pos */
        /* get initial lnprob values */
        iwalker=0;
        for(rank=1; rank<nprocs; rank++){
            current_task[rank]=iwalker;
            MPI_Send(&start_pos[iwalker].pars[0], npars, MPI_DOUBLE, rank, WORKTAG, MPI_COMM_WORLD);
            iwalker++;
        }

        /* receive lnprob's as they come in, and send out remaining walker positions */
        while(iwalker<nwalkers){
#ifndef WRITE_EXTRA_DOUBLES
            MPI_Recv(&lnprob_tmp, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
#else
            MPI_Recv(&extra_doubles[0], nextra+1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            lnprob_tmp=extra_doubles[nextra];
#endif
            irecv=current_task[status.MPI_SOURCE];
            start_pos[irecv].lnprob = lnprob_tmp;
            start_pos[irecv].rank = status.MPI_SOURCE;
            current_task[status.MPI_SOURCE] = iwalker;
#ifdef WRITE_EXTRA_DOUBLES
            for(iextra=0; iextra<nextra; iextra++){
                start_pos[irecv].extra_doubles[iextra]=extra_doubles[iextra];
            }
#endif
            MPI_Send(&start_pos[iwalker].pars[0], npars, MPI_DOUBLE, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
            iwalker++;
        }

        /* receive remaining lnprobs -- should be one for each worker, but they arrive in unknown order */
        for(rank = 1; rank < nprocs; rank++) {
#ifndef WRITE_EXTRA_DOUBLES
            MPI_Recv(&lnprob_tmp, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
#else
            MPI_Recv(&extra_doubles[0], nextra+1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            lnprob_tmp=extra_doubles[nextra];
#endif
            irecv=current_task[status.MPI_SOURCE];
            start_pos[irecv].rank = status.MPI_SOURCE;
            start_pos[irecv].lnprob = lnprob_tmp;
#ifdef WRITE_EXTRA_DOUBLES
            for(iextra=0; iextra<nextra; iextra++){
                start_pos[irecv].extra_doubles[iextra]=extra_doubles[iextra];
            }
#endif
        }

        /* fill the two ensembles with starting positions */
        for(iwalker=0; iwalker<nwalkers_over_two; iwalker++){
            ensemble_A[iwalker].accept=1;
            ensemble_B[iwalker].accept=1;
            ensemble_A[iwalker].lnprob=start_pos[iwalker].lnprob;
            ensemble_B[iwalker].lnprob=start_pos[iwalker+nwalkers_over_two].lnprob;
            ensemble_A[iwalker].rank=start_pos[iwalker].rank;
            ensemble_B[iwalker].rank=start_pos[iwalker+nwalkers_over_two].rank;
            for(ipar=0; ipar<npars; ipar++){
                ensemble_A[iwalker].pars[ipar]=start_pos[iwalker].pars[ipar];
                ensemble_B[iwalker].pars[ipar]=start_pos[iwalker+nwalkers_over_two].pars[ipar];
            }
#ifdef WRITE_EXTRA_DOUBLES
            for(iextra=0; iextra<nextra; iextra++){
                ensemble_A[iwalker].extra_doubles[iextra]=start_pos[iwalker].extra_doubles[iextra];
                ensemble_B[iwalker].extra_doubles[iextra]=start_pos[iwalker+nwalkers_over_two].extra_doubles[iextra];
            }
#endif
        }

        write_header(nsteps, nwalkers, npars, nextra, fname);
        write_step(npars, nwalkers_over_two, 0, nextra, ensemble_A, ensemble_B, fname);
        istart = 1;
    }
    else{
        /* read the last nwalkers lines of the output file in order to get start_pos */
        Nlines = getNlines(fname,'%'); // using the wrong comment on purpose to count all lines
        startline = Nlines - nwalkers;
        iline = 0;

        file = fopen(fname,"r");
        while(iline<Nlines){
            /* I should implement a smarter way to just skip over the lines */
            if(iline<startline){
                fgets(buffer, MAXLINESIZE,file);
            }
            else if(iline<(startline+nwalkers_over_two)){
                iwalker = iline - startline;
                fscanf(file,"%*d %*d %d %lf", &ensemble_A[iwalker].accept, &ensemble_A[iwalker].lnprob);
                for(ipar=0; ipar<npars; ipar++){
                    fscanf(file,"%lf", &ensemble_A[iwalker].pars[ipar]);
                }
#ifdef WRITE_EXTRA_DOUBLES
                for(iextra=0; iextra<nextra; iextra++){
                    fscanf(file,"%lf", &ensemble_A[iwalker].extra_doubles[iextra]);
                }
#endif
            }
            else{
                iwalker = iline - startline - nwalkers_over_two;
                fscanf(file,"%*d %*d %d %lf", &ensemble_B[iwalker].accept, &ensemble_B[iwalker].lnprob);
                for(ipar=0; ipar<npars; ipar++){
                    fscanf(file,"%lf", &ensemble_B[iwalker].pars[ipar]);
                }
#ifdef WRITE_EXTRA_DOUBLES
                for(iextra=0; iextra<nextra; iextra++){
                    fscanf(file,"%lf", &ensemble_B[iwalker].extra_doubles[iextra]);
                }
#endif
            }
            iline++;
        }
        fclose(file);

        offset = Nlines/nwalkers; //Nlines is too big but OK bc int division
        istart = offset;
        nsteps += offset;
    }

    /* begin the chain */
    for(istep=istart; istep<nsteps; istep++){
        /*----------------------------- ensemble A ---------------------------*/
        create_trials(nwalkers_over_two, npars, a, z_array, trial, ensemble_A, ensemble_B);
        /* send out first batch of trial positions */
        iwalker=0;
        for(rank=1; rank<nprocs; rank++){
            current_task[rank]=iwalker;
            MPI_Send(&trial[iwalker].pars[0], npars, MPI_DOUBLE, rank, WORKTAG, MPI_COMM_WORLD);
            iwalker++;
        }

        /* receive lnprob's as they come in, and send out remaining walker positions */
        while(iwalker<nwalkers_over_two){
#ifndef WRITE_EXTRA_DOUBLES
            MPI_Recv(&lnprob_tmp, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
#else
            MPI_Recv(&extra_doubles[0], nextra+1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            lnprob_tmp=extra_doubles[nextra];
#endif
            irecv = current_task[status.MPI_SOURCE];
            trial[irecv].lnprob = lnprob_tmp;
            ensemble_A[irecv].rank = status.MPI_SOURCE;
            current_task[status.MPI_SOURCE] = iwalker;
#ifdef WRITE_EXTRA_DOUBLES
            for(iextra=0; iextra<nextra; iextra++){
                trial[irecv].extra_doubles[iextra]=extra_doubles[iextra];
            }
#endif
            MPI_Send(&trial[iwalker].pars[0], npars, MPI_DOUBLE, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
            iwalker++;
        }

        /* receive remaining lnprobs -- should be one for each worker, but they arrive in unknown order */
        for (rank = 1; rank < nprocs; rank++) {
#ifndef WRITE_EXTRA_DOUBLES
            MPI_Recv(&lnprob_tmp, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
#else
            MPI_Recv(&extra_doubles[0], nextra+1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            lnprob_tmp=extra_doubles[nextra];
#endif
            irecv = current_task[status.MPI_SOURCE];
            trial[irecv].lnprob = lnprob_tmp;
            ensemble_A[irecv].rank = status.MPI_SOURCE;
#ifdef WRITE_EXTRA_DOUBLES
            for(iextra=0; iextra<nextra; iextra++){
                trial[irecv].extra_doubles[iextra]=extra_doubles[iextra];
            }
#endif
        }

        /* update walker positions */
        update_positions(nwalkers_over_two, npars, nextra, z_array, ensemble_A, trial);

        /*----------------------------- ensemble B ---------------------------*/
        create_trials(nwalkers_over_two, npars, a, z_array, trial, ensemble_B, ensemble_A);

        /* send out first batch of trial positions */
        iwalker=0;
        for(rank=1; rank<nprocs; rank++){
            current_task[rank] = iwalker;
            MPI_Send(&trial[iwalker].pars[0], npars, MPI_DOUBLE, rank, WORKTAG, MPI_COMM_WORLD);
            iwalker++;
        }

        /* receive lnprob's as they come in, and send out remaining walker positions */
        while(iwalker<nwalkers_over_two){
#ifndef WRITE_EXTRA_DOUBLES
            MPI_Recv(&lnprob_tmp, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
#else
            MPI_Recv(&extra_doubles[0], nextra+1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            lnprob_tmp=extra_doubles[nextra];
#endif
            irecv = current_task[status.MPI_SOURCE];
            trial[irecv].lnprob = lnprob_tmp;
            ensemble_B[irecv].rank = status.MPI_SOURCE;
            current_task[status.MPI_SOURCE] = iwalker;
#ifdef WRITE_EXTRA_DOUBLES
            for(iextra=0; iextra<nextra; iextra++){
                trial[irecv].extra_doubles[iextra]=extra_doubles[iextra];
            }
#endif
            MPI_Send(&trial[iwalker].pars[0], npars, MPI_DOUBLE, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
            iwalker++;
        }

        /* receive remaining lnprobs -- should be one for each worker, but they arrive in unknown order */
        for (rank = 1; rank < nprocs; rank++) {
#ifndef WRITE_EXTRA_DOUBLES
            MPI_Recv(&lnprob_tmp, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
#else
            MPI_Recv(&extra_doubles[0], nextra+1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            lnprob_tmp=extra_doubles[nextra];
#endif
            irecv = current_task[status.MPI_SOURCE];
            trial[irecv].lnprob = lnprob_tmp;
            ensemble_B[irecv].rank = status.MPI_SOURCE;
#ifdef WRITE_EXTRA_DOUBLES
            for(iextra=0; iextra<nextra; iextra++){
                trial[irecv].extra_doubles[iextra]=extra_doubles[iextra];
            }
#endif
        }

        /* update walker positions */
        update_positions(nwalkers_over_two, npars, nextra, z_array, ensemble_B, trial);
        write_step(npars, nwalkers_over_two, istep, nextra, ensemble_A, ensemble_B, fname);
    }

    /* chain complete -- tell workers to exit */
    for (rank = 1; rank < nprocs; rank++) {
        MPI_Send(0, 0, MPI_DOUBLE, rank, DIETAG, MPI_COMM_WORLD);
    }
    free(z_array);
    free(current_task);
#ifdef WRITE_EXTRA_DOUBLES
    free(extra_doubles);
#endif
    free_walkers(nwalkers_over_two, ensemble_A);
    free_walkers(nwalkers_over_two, ensemble_B);
    free_walkers(nwalkers_over_two, trial);
}

/* -------------------------------------------------------------------------- */
#ifndef WRITE_EXTRA_DOUBLES
void worker(int npars, const void *userdata, double (*lnprob)(const double *, int, const void *))
{
    double lnprob_new;
    double *pars;
    MPI_Status status;

    pars=calloc(npars, sizeof(double));
    while(1){
        MPI_Recv(&pars[0], npars, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == DIETAG) {
            break;
        }
        lnprob_new = lnprob(pars, npars, userdata);
        MPI_Send(&lnprob_new, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    free(pars);
}
#else
void worker(int npars, int nextra, const void *userdata,
            double (*lnprob)(const double *, int, const void *, int, double *))
{
    double lnprob_new;
    double *pars,*extra_doubles;
    MPI_Status status;

    extra_doubles=calloc(nextra+1,sizeof(double));
    pars=calloc(npars, sizeof(double));
    while(1){
        MPI_Recv(&pars[0], npars, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == DIETAG) {
            break;
        }
        lnprob_new = lnprob(pars, npars, userdata, nextra, extra_doubles);
        extra_doubles[nextra]=lnprob_new;
        MPI_Send(&extra_doubles[0], nextra+1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    free(pars);
    free(extra_doubles);
}
#endif

/* -------------------------------------------------------------------------- */
#ifndef WRITE_EXTRA_DOUBLES
void run_chain(int *argc, char ***argv, int nwalkers, int nsteps, int npars,
               int nburn, int resume, double a, walker_pos *start_pos,
               double (*lnprob)(const double *, int, const void *),
               const void *userdata, const char *fname)
{
    // first check that the file has more lines than walkers if we are resuming a chain
    size_t nlines;
    if(resume==1){
        if(access(fname,R_OK)==-1){
            fprintf(stderr, "Cannot resume chain from file %s: no read access\n",
                    fname);
            return;
        }
        nlines = getNlines(fname, '#');
        if(nlines<nwalkers){
            fprintf(stderr, "Cannot resume chain from file %s: has %zu lines, less than nwalkers\n",
                    fname, nlines);
            return;
        }
    }

    /*========================================================================*/
    /* MPI stuff */
    int nprocs, rank, nworkers;
    int user_control, flag;
    user_control = 0;
    MPI_Initialized(&flag);

    if (flag){
        user_control = 1; //the user has already called MPI_Init; will be responsible for calling MPI_Finalize
    }
    else{
        MPI_Init(argc,argv);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* chain things */
    int nwalkers_over_two;

    /*========================================================================*/
    /* check that nwalkers is even */
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
            fprintf(stderr, "Attempting to split a 'half-ensemble' of %d walkers across %d workers.\n",
                    nwalkers_over_two, nworkers);
            fprintf(stderr, "Change the number of workers and walkers\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    if(rank==0){
        manager(nwalkers, nsteps, npars, 0, nburn, resume, a, start_pos, fname);
    }
    else worker(npars, userdata, lnprob);

    /* end MPI */
    if (user_control==0) MPI_Finalize();
}

#else
void run_chain(int *argc, char ***argv, int nwalkers, int nsteps, int npars,
               int nextra, int nburn, int resume, double a, walker_pos *start_pos,
               double (*lnprob)(const double *, int, const void *, int, double *),
               const void *userdata, const char *fname)
{
    // first check that the file has more lines than walkers if we are resuming a chain
    size_t nlines;
    if(resume==1){
        if(access(fname,R_OK)==-1){
            fprintf(stderr, "Cannot resume chain from file %s: no read access\n",
                    fname);
            return;
        }
        nlines = getNlines(fname, '#');
        if(nlines<nwalkers){
            fprintf(stderr, "Cannot resume chain from file %s: has %zu lines, less than nwalkers\n",
                    fname, nlines);
            return;
        }
    }
    if(nextra<1){
        fprintf(stderr, "WRITE_EXTRA_DOUBLES enabled but nextra=%d. Must be > 0\n", nextra);
        fprintf(stderr, "Disable WRITE_EXTRA_DOUBLES in Makefile if you have no extra doubles to write\n");
        exit(EXIT_FAILURE);
    }

    /*========================================================================*/
    /* MPI stuff */
    int nprocs, rank, nworkers;
    int user_control, flag;
    user_control = 0;
    MPI_Initialized(&flag);

    if (flag){
        user_control = 1; //the user has already called MPI_Init; will be responsible for calling MPI_Finalize
    }
    else{
        MPI_Init(argc,argv);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* chain things */
    int nwalkers_over_two;

    /*========================================================================*/
    /* check that nwalkers is even */
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
            fprintf(stderr, "Attempting to split a 'half-ensemble' of %d walkers across %d workers.\n",
                    nwalkers_over_two, nworkers);
            fprintf(stderr, "Change the number of workers and walkers\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    if(rank==0){
        manager(nwalkers, nsteps, npars, nextra, nburn, resume, a, start_pos, fname);
    }
    else worker(npars, nextra, userdata, lnprob);

    /* end MPI */
    if (user_control==0) MPI_Finalize();
}
#endif