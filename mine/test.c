#include "emcee.h"

int main( int argc, char ** argv )
{
    /*========================================================================*/
    /* MPI stuff */
    int nprocs, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* generic indices */
    int i, j;

    /* chain things */
    int nsteps, nwalkers, npars, nwalkers_over_two;
    ensemble *my_ensemble;
    walker_pos *my_walkers;
    chain *my_chain;

    /* more MPI stuff */
    int slice_length, lower_ind, upper_ind, remain;
    int mpi_disp[nprocs], counts[nprocs];
    int tmp_start, tmp_slice;

    /*========================================================================*/

    /* Set up basic chain parameters */
    nsteps=(int)NSTEPS;
    nwalkers=(int)NWALKERS;
    npars=(int)NPARS;
    nwalkers_over_two=nwalkers/2;

    if(nprocs>nwalkers_over_two){
        if(rank==0){
            fprintf(stderr, "Attempting to split a 'half-ensemble' of %d walkers across %d processes is silly.\n",
                    nwalkers_over_two, nprocs);
            fprintf(stderr, "I don't know how to do that. Change the values in the 'pars.h' file\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
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
    remain=nwalkers%nprocs;

    /* Fill mpi_disp and counts for MPI comm */
    for(i=0; i<nprocs; i++)
    {
        slice_length=nwalkers/nprocs;
        lower_ind=i*slice_length;
        if(i<remain){
            lower_ind+=i;
            slice_length++;
        }
        else{
            lower_ind+=remain;
        }
        mpi_disp[i]=lower_ind;
        counts[i]=slice_length;
    }

    /* Have each rank set its slice indices */
    lower_ind = mpi_disp[rank];
    slice_length = counts[rank];
    upper_ind = lower_ind + slice_length;

    /*========================================================================*/

    /* set up chain, ensemble, walkers */
    if (rank==0) my_chain=allocate_chain(nsteps,nwalkers,npars);

    my_ensemble=allocate_ensemble(nwalkers,npars);
    my_walkers=allocate_walkers(slice_length);

    /*========================================================================*/

    /* fill with some data */
    for(i=0; i<slice_length; i++){
        my_walkers[i].accept = 15;
        my_walkers[i].lnprob = 279.6;
        for(j=0; j<npars; j++){
            my_walkers[i].pars[j] = (double)(j)+342.1;
        }
    }

    /* gather data from each proc */
    MPI_Allgatherv(&my_walkers[0], slice_length, MPI_WALKER,
                   &my_ensemble->walker[0], counts, mpi_disp,
                   MPI_WALKER, MPI_COMM_WORLD);

    /* check allgatherv results */
    for(i=0; i<nwalkers; i++){
        if(my_ensemble->walker[i].accept != 15){
            fprintf(stderr, "Error in accept: rank %d walker %d\n",rank,i);
        }
        if(my_ensemble->walker[i].lnprob != 279.6){
            fprintf(stderr, "Error in lnprob: rank %d walker %d\n",rank,i);
        }
        for(j=0; j<npars; j++){
            if(my_ensemble->walker[i].pars[j] != (double)(j)+342.1){
                fprintf(stderr, "Error in pars: rank %d walker %d par %d\n",rank,i,j);
            }
        }
    }

    /*========================================================================*/

    /* free chain */
    if(rank==0) free_chain(my_chain);
    free_ensemble(my_ensemble);
    free_walkers(my_walkers);

    /* end MPI */
    MPI_Type_free(&MPI_WALKER);
    MPI_Finalize();
    return 0;
}