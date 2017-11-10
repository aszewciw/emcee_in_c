#include "emcee.h"

// struct mydata {
//     size_t ndata;
//     const double *data;
//     double ivar; // same error for each
// };

// double lnprob(const double *pars, size_t npars, const void *userdata)
// {
//     double chi2=0, diff=0;
//     const struct mydata *mydata = userdata;

//     for (size_t i=0; i<mydata->ndata; i++) {
//         diff = mydata->data[i]-pars[0];
//         chi2 += diff*diff*mydata->ivar;
//     }

//     double lnprob = -0.5*chi2;

//     return lnprob;
// }

void run_chain(int *argc, char ***argv, double *centers, double *widths){
    /*========================================================================*/
    /* MPI stuff */
    int nprocs, rank;
    MPI_Init(argc,argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* generic indices */
    int iproc, iwalker, ipar, istep;

    /* chain things */
    int nsteps, nwalkers, npars, nwalkers_over_two;
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

    /* set up the rng here */
    srand((unsigned) (time(&t)+rank));
    fprintf(stderr, "process: %d, rng: %d\n", rand());

    /* Have each process make its own guess. We'll overwrite it in a second */
    start_pos = make_guess(centers,widths,nwalkers,npars);

    int current_rank=0;
    while(current_rank)<nprocs{
        if(rank==current_rank){
            fprintf(stderr, "Rank %d\n", rank);
            for(iwalker=0;iwalker<nwalkers;iwalker++){
                for(ipar=0;ipar<npars;ipar++){
                    fprintf(stderr, "\t%lf\n", start_pos[iwalker].pars[ipar]);
                }
            }
        }
        current_rank++;
    }

    /* Have process 0 send data to all others */
    MPI_Bcast(&start_pos[0], nwalkers, MPI_WALKER, 0, MPI_COMM_WORLD);
    if(rank==0) fprintf(stderr, "\n");
    while(current_rank)<nprocs{
        if(rank==current_rank){
            fprintf(stderr, "Rank %d\n", rank);
            for(iwalker=0;iwalker<nwalkers;iwalker++){
                for(ipar=0;ipar<npars;ipar++){
                    fprintf(stderr, "\t%lf\n", start_pos[iwalker].pars[ipar]);
                }
            }
        }
        current_rank++;
    }

    /* fill with some data */
    // for(iwalker=0; iwalker<slice_length; iwalker++){
    //     my_walkers[iwalker].accept = 15;
    //     my_walkers[iwalker].lnprob = 279.6;
    //     for(ipar=0; ipar<npars; ipar++){
    //         my_walkers[iwalker].pars[ipar] = (double)(ipar)+342.1;
    //     }
    // }

    /* gather data from each proc */
    // MPI_Allgatherv(&my_walkers[0], slice_length, MPI_WALKER,
    //                &ensemble_A->walker[0], counts, mpi_disp,
    //                MPI_WALKER, MPI_COMM_WORLD);

    // /* check allgatherv results */
    // for(i=0; i<nwalkers; i++){
    //     if(my_ensemble->walker[i].accept != 15){
    //         fprintf(stderr, "Error in accept: rank %d walker %d\n",rank,i);
    //     }
    //     if(my_ensemble->walker[i].lnprob != 279.6){
    //         fprintf(stderr, "Error in lnprob: rank %d walker %d\n",rank,i);
    //     }
    //     for(j=0; j<npars; j++){
    //         if(my_ensemble->walker[i].pars[j] != (double)(j)+342.1){
    //             fprintf(stderr, "Error in pars: rank %d walker %d par %d\n",rank,i,j);
    //         }
    //     }
    // }

    /*========================================================================*/

    /* free chain */
    if(rank==0) free_chain(my_chain);
    free_ensemble(my_ensemble);
    free_walkers(my_walkers);

    /* end MPI */
    MPI_Type_free(&MPI_WALKER);
    MPI_Finalize();
}

int main( int argc, char ** argv )
{
    double centers[NPARS] = {1.0,2.0};
    double widths[NPARS] = {1.0,1.0}

    run_chain(&argc, &argv, centers, widths);
    return 0;
}















// int main( int argc, char ** argv )
// {
//     /*========================================================================*/
//     /* MPI stuff */
//     int nprocs, rank;
//     MPI_Init(&argc,&argv);
//     MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//     /* generic indices */
//     int i, j;

//     /* chain things */
//     int nsteps, nwalkers, npars, nwalkers_over_two;
//     ensemble *my_ensemble;
//     walker_pos *my_walkers;
//     chain *my_chain;

//     /* more MPI stuff */
//     int slice_length, lower_ind, upper_ind, remain;
//     int mpi_disp[nprocs], counts[nprocs];
//     int tmp_start, tmp_slice;

//     /*========================================================================*/

//     /* Set up basic chain parameters */
//     nsteps=(int)NSTEPS;
//     nwalkers=(int)NWALKERS;
//     npars=(int)NPARS;
//     nwalkers_over_two=nwalkers/2;

//     if(nprocs>nwalkers_over_two){
//         if(rank==0){
//             fprintf(stderr, "Attempting to split a 'half-ensemble' of %d walkers across %d processes is silly.\n",
//                     nwalkers_over_two, nprocs);
//             fprintf(stderr, "I don't know how to do that. Change the values in the 'pars.h' file\n");
//         }
//         MPI_Barrier(MPI_COMM_WORLD);
//         MPI_Finalize();
//         exit(EXIT_FAILURE);
//     }

//     /*========================================================================*/

//     /* make an MPI custom structure */
//     int blocklen[3] = {1, 1, npars};
//     MPI_Datatype MPI_WALKER;
//     MPI_Datatype type[3] = { MPI_INT, MPI_DOUBLE, MPI_DOUBLE };
//     MPI_Aint disp[3];

//     disp[0] = offsetof(walker_pos,accept);
//     disp[1] = offsetof(walker_pos,lnprob);
//     disp[2] = offsetof(walker_pos,pars);
//     MPI_Type_create_struct(3,blocklen,disp,type,&MPI_WALKER);
//     MPI_Type_commit(&MPI_WALKER);

//     /*========================================================================*/

//     /* Establish slice of walkers for each process to handle */
//     remain=nwalkers%nprocs;

//     /* Fill mpi_disp and counts for MPI comm */
//     for(i=0; i<nprocs; i++)
//     {
//         slice_length=nwalkers/nprocs;
//         lower_ind=i*slice_length;
//         if(i<remain){
//             lower_ind+=i;
//             slice_length++;
//         }
//         else{
//             lower_ind+=remain;
//         }
//         mpi_disp[i]=lower_ind;
//         counts[i]=slice_length;
//     }

//     /* Have each rank set its slice indices */
//     lower_ind = mpi_disp[rank];
//     slice_length = counts[rank];
//     upper_ind = lower_ind + slice_length;

//     /*========================================================================*/

//     /* set up chain, ensemble, walkers */
//     if (rank==0) my_chain=allocate_chain(nsteps,nwalkers,npars);

//     my_ensemble=allocate_ensemble(nwalkers,npars);
//     my_walkers=allocate_walkers(slice_length);

//     /*========================================================================*/

//     /* fill with some data */
//     for(i=0; i<slice_length; i++){
//         my_walkers[i].accept = 15;
//         my_walkers[i].lnprob = 279.6;
//         for(j=0; j<npars; j++){
//             my_walkers[i].pars[j] = (double)(j)+342.1;
//         }
//     }

//     /* gather data from each proc */
//     MPI_Allgatherv(&my_walkers[0], slice_length, MPI_WALKER,
//                    &my_ensemble->walker[0], counts, mpi_disp,
//                    MPI_WALKER, MPI_COMM_WORLD);

//     /* check allgatherv results */
//     for(i=0; i<nwalkers; i++){
//         if(my_ensemble->walker[i].accept != 15){
//             fprintf(stderr, "Error in accept: rank %d walker %d\n",rank,i);
//         }
//         if(my_ensemble->walker[i].lnprob != 279.6){
//             fprintf(stderr, "Error in lnprob: rank %d walker %d\n",rank,i);
//         }
//         for(j=0; j<npars; j++){
//             if(my_ensemble->walker[i].pars[j] != (double)(j)+342.1){
//                 fprintf(stderr, "Error in pars: rank %d walker %d par %d\n",rank,i,j);
//             }
//         }
//     }

//     /*========================================================================*/

//     /* free chain */
//     if(rank==0) free_chain(my_chain);
//     free_ensemble(my_ensemble);
//     free_walkers(my_walkers);

//     /* end MPI */
//     MPI_Type_free(&MPI_WALKER);
//     MPI_Finalize();
//     return 0;
// }