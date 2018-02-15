#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <float.h>
#include <math.h>
#include <string.h>

#define WORKTAG     1
#define DIETAG     2

void master(int ntasks)
{

    int nprocs, rank, work, itask, tmpres, ires;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int results[ntasks];
    for(itask=0;itask<ntasks;itask++){
        results[itask]=0;
    }
    itask=0;
    for (rank = 1; rank < nprocs; rank++) {
        work = itask;
        MPI_Send(&work,1,MPI_INT,rank,WORKTAG,MPI_COMM_WORLD);
        itask++;
    }

    ires=0;
    while (itask<ntasks) {
        work = itask;
        MPI_Recv(&tmpres, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        results[ires]=tmpres;
        ires++;
        MPI_Send(&work, 1, MPI_INT, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
        itask++;
    }

/*
* Receive results for outstanding work requests.
*/
    for (rank = 1; rank < nprocs; rank++) {
        MPI_Recv(&tmpres, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        results[ires]=tmpres;
        ires++;
    }

/*
* Tell all the slaves to exit.
*/
    for(ires=0;ires<ntasks;ires++){
        fprintf(stderr, "\t\tRank 0: %d\n", results[ires]);
    }

    for (rank = 1; rank < nprocs; rank++) {
        MPI_Send(0, 1, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
    }

}

void slave(int myrank)
{
    int work,result;
    MPI_Status status;
    while(1) {
        MPI_Recv(&work, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        fprintf(stderr, "rank %d; work %d\n", myrank,work);
        if (status.MPI_TAG == DIETAG) {
            break;
        }
        result = myrank;
        MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    return;
}

int main(int argc, char **argv)
{
    int myrank,ntasks;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    if(argc!=2){
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    ntasks=atoi(argv[1]);

    if (myrank == 0) {
        master(ntasks);
    }
    else {
        slave(myrank);
    }
    MPI_Finalize();
}