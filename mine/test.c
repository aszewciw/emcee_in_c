#include "emcee.h"

struct mydata {
    int ndata;
    const double *data;
    double ivar; // same error for each
};

double lnprob(const double *pars, int npars, const void *userdata)
{
    int i;
    double chi2,diff,lnprob;
    const struct mydata *mydata = userdata;

    chi2=0;
    diff=0;

    for (i=0; i<mydata->ndata; i++) {
        diff = mydata->data[i]-pars[0];
        chi2 += diff*diff*mydata->ivar;
    }

    lnprob = -0.5*chi2;

    return lnprob;
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

    // int current_rank=0;

    // while(current_rank<nprocs){
    //     if (rank==current_rank){
    //         fprintf(stderr, "\n\nRank: %d\n", rank);
    //         for(iwalker=0;iwalker<nwalkers_over_two;iwalker++){
    //             fprintf(stderr, "Walker %d, accept %d, lnprob %lf\n",
    //                     iwalker, ensemble_B->walker[iwalker].accept,
    //                     ensemble_B->walker[iwalker].lnprob);
    //             for(ipar=0;ipar<npars;ipar++){
    //                 fprintf(stderr, "\t%lf", ensemble_B->walker[iwalker].pars[ipar]);
    //             }
    //             fprintf(stderr, "\n");
    //         }
    //     }
    //     current_rank++;
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    /* gather data from each proc */
    MPI_Allgatherv(&my_walkers[0], slice_length, MPI_WALKER,
                   &ensemble_B->walker[0], counts, mpi_disp,
                   MPI_WALKER, MPI_COMM_WORLD);

    // current_rank=0;
    // while(current_rank<nprocs){
    //     if (rank==current_rank){
    //         fprintf(stderr, "\n\nRank: %d\n", rank);
    //         for(iwalker=0;iwalker<nwalkers_over_two;iwalker++){
    //             fprintf(stderr, "Walker %d, accept %d, lnprob %lf\n",
    //                     iwalker, ensemble_B->walker[iwalker].accept,
    //                     ensemble_B->walker[iwalker].lnprob);
    //             for(ipar=0;ipar<npars;ipar++){
    //                 fprintf(stderr, "\t%lf", ensemble_B->walker[iwalker].pars[ipar]);
    //             }
    //             fprintf(stderr, "\n");
    //         }
    //     }
    //     current_rank++;
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

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

    // current_rank=0;
    // while(current_rank<nprocs){
    //     if (rank==current_rank){
    //         fprintf(stderr, "\n\nRank: %d\n", rank);
    //         for(iwalker=0;iwalker<nwalkers_over_two;iwalker++){
    //             fprintf(stderr, "Walker %d, accept %d, lnprob %lf\n",
    //                     iwalker, ensemble_A->walker[iwalker].accept,
    //                     ensemble_A->walker[iwalker].lnprob);
    //             for(ipar=0;ipar<npars;ipar++){
    //                 fprintf(stderr, "\t%lf", ensemble_A->walker[iwalker].pars[ipar]);
    //             }
    //             fprintf(stderr, "\n");
    //         }
    //     }
    //     current_rank++;
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    /* gather data from each proc */
    MPI_Allgatherv(&my_walkers[0], slice_length, MPI_WALKER,
                   &ensemble_A->walker[0], counts, mpi_disp,
                   MPI_WALKER, MPI_COMM_WORLD);

    // current_rank=0;
    // while(current_rank<nprocs){
    //     if (rank==current_rank){
    //         fprintf(stderr, "\n\nRank: %d\n", rank);
    //         for(iwalker=0;iwalker<nwalkers_over_two;iwalker++){
    //             fprintf(stderr, "Walker %d, accept %d, lnprob %lf\n",
    //                     iwalker, ensemble_A->walker[iwalker].accept,
    //                     ensemble_A->walker[iwalker].lnprob);
    //             for(ipar=0;ipar<npars;ipar++){
    //                 fprintf(stderr, "\t%lf", ensemble_A->walker[iwalker].pars[ipar]);
    //             }
    //             fprintf(stderr, "\n");
    //         }
    //     }
    //     current_rank++;
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

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

int main( int argc, char ** argv )
{
    int ndata=10;
    int i;
    double a=2.0;

    double truepars[1] = {1};
    double guess[1] = {0};
    double ballsize[1] = {0};
    double fracerr=0.1;

    double err=fracerr*truepars[0];

    double *data = malloc(ndata*sizeof(double));
    for(i=0; i<ndata; i++) {
        data[i] = truepars[0] + err*mca_randn();
    }

    struct mydata mydata;
    mydata.ndata = ndata;
    mydata.data = (const double*) data;
    mydata.ivar = 1/(err*err);

    guess[0] = truepars[0] + err*mca_randn();
    ballsize[0] = 1.0;

    const char fname[]="chain.dat";
    run_chain(&argc, &argv, a, guess, ballsize, &lnprob, &mydata, fname);
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