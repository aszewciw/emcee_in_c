#include "emcee.h"

int main( int argc, char ** argv )
{
    int nprocs,rank;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int i,j;
    int nsteps=(int)NSTEPS;
    int nwalkers=(int)NWALKERS;
    int npars=(int)NPARS;
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
    walker_pos *my_walkers=allocate_walkers(slice_length);

    // have each chain fill its walker_pos with nonsense that we'll check later
    for(i=0; i<slice_length; i++){
        my_walkers[i].accept = 15;
        my_walkers[i].lnprob = 279.6;
        for(j=0; j<npars; j++){
            my_walkers[i].pars[j] = (double)(j)+342.1;
        }
    }

    int current_rank=0;

    MPI_Datatype MPI_WALKER;
    MPI_Datatype type[3] = { MPI_INT, MPI_DOUBLE, MPI_DOUBLE };
    int blocklen[3] = { 1, 1, npars }; // size of each data element in struct

    MPI_Aint disp[3]; // array of displacements; one for each data member
    disp[0] = offsetof(walker_pos,accept);
    disp[1] = offsetof(walker_pos,lnprob);
    disp[2] = offsetof(walker_pos,pars);

    MPI_Type_create_struct(3,blocklen,disp,type,&MPI_WALKER);
    MPI_Type_commit(&MPI_WALKER);

    int mpi_disp[nprocs];
    int counts[nprocs];
    int tmp_start, tmp_slice;

    for(i=0; i<nprocs; i++)
    {
        tmp_slice=nwalkers/nprocs;
        tmp_start=i*slice_length;
        if(i<remain){
            tmp_start+=i;
            tmp_slice++;
        }
        else{
            tmp_start+=remain;
        }
        mpi_disp[i]=tmp_start;
        counts[i]=tmp_slice;
    }


    MPI_Allgatherv(&my_walkers[0], slice_length, MPI_WALKER,
                 &my_ensemble->walker[0], counts, mpi_disp,
                 MPI_WALKER, MPI_COMM_WORLD);

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

    if(rank==0) free_chain(my_chain);
    free_ensemble(my_ensemble);
    free_walkers(my_walkers);

    MPI_Type_free(&MPI_WALKER);
    MPI_Finalize();
    return 0;
}