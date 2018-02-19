#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <float.h>
#include <math.h>
#include <string.h>

#define WORKTAG 1
#define DIETAG 2
#define NPARS 3

typedef struct walker_pos{
    int accept;
    double lnprob;
    double pars[NPARS];
} walker_pos;

walker_pos* allocate_walkers(size_t nwalkers){
    struct walker_pos *self=calloc(nwalkers,sizeof(walker_pos));
    if (self==NULL) {
        fprintf(stderr,"Could not allocate struct walker_pos\n");
        exit(EXIT_FAILURE);
    }
    return self;
}

void manager(int ntasks, MPI_Datatype MPI_WALKER)
{

    int nprocs, rank, work, itask, tmpres, ires;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    walker_pos *my_walkers;

    my_walkers = allocate_walkers(1);

    // double results[ntasks];
    double lnprob[ntasks];
    int current_task[nprocs];     // array of each proc's currently assigned task
    // for(itask=0;itask<ntasks;itask++){
    //     results[itask]=0;
    // }
    for(itask=0;itask<ntasks;itask++){
        lnprob[itask]=0;
    }
    for(rank=0; rank<nprocs; rank++){
        current_task[rank]=-1;
    }
    itask=0;
    for (rank = 1; rank < nprocs; rank++) {
        current_task[rank]=itask;
        // MPI_Send(&work,1,MPI_INT,rank,WORKTAG,MPI_COMM_WORLD);
        MPI_Send(&my_walkers[0],1,MPI_WALKER,rank,WORKTAG,MPI_COMM_WORLD);
        itask++;
    }

    while (itask<ntasks) {
        // work = itask;
        // MPI_Recv(&tmpres, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        // results[ires]=tmpres;
        MPI_Recv(&my_walkers[0], 1, MPI_WALKER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        ires=current_task[status.MPI_SOURCE];
        lnprob[ires]=my_walkers[0].lnprob;
        current_task[status.MPI_SOURCE]=itask;
        // MPI_Send(&work, 1, MPI_INT, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
        MPI_Send(&my_walkers[0], 1, MPI_WALKER, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
        itask++;
    }

/*
* Receive results for outstanding work requests.
*/
    for (rank = 1; rank < nprocs; rank++) {
        // MPI_Recv(&tmpres, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        // results[ires]=tmpres;
        MPI_Recv(&my_walkers[0], 1, MPI_WALKER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        ires=current_task[status.MPI_SOURCE];
        lnprob[ires]=my_walkers[0].lnprob;
    }

/*
* Tell all the workers to exit.
*/
    for(ires=0;ires<ntasks;ires++){
        // fprintf(stderr, "\t\tRank 0: %d\n", results[ires]);
        fprintf(stderr, "lnprob %d: %lf\n", ires, lnprob[ires]);
    }
    for(rank=0; rank<nprocs; rank++){
        fprintf(stderr, "curr task: %d\n", current_task[rank]);
    }
    for (rank = 1; rank < nprocs; rank++) {
        // MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
        MPI_Send(0, 0, MPI_WALKER, rank, DIETAG, MPI_COMM_WORLD);
    }

    free(my_walkers);

}

void worker(int myrank, MPI_Datatype MPI_WALKER)
{
    // int work,result;
    double result;
    int ipar,i,max_ind;
    walker_pos *my_walkers, *manager_walkers;
    my_walkers = allocate_walkers(1);
    manager_walkers = allocate_walkers(1);

    if(myrank==3){
        max_ind=1000000000;
    }
    else{
        max_ind=10;
    }
    fprintf(stderr, "%d\n", max_ind);
    my_walkers[0].accept=myrank;
    my_walkers[0].lnprob=(double)myrank;
    for(ipar=0;ipar<NPARS;ipar++){
        my_walkers[0].pars[ipar]=myrank;
    }
    MPI_Status status;
    int work_counter=0;
    while(1) {
        // MPI_Recv(&work, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        // fprintf(stderr, "rank %d; work %d\n", myrank,work);
        MPI_Recv(&manager_walkers[0], 1, MPI_WALKER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == DIETAG) {
            break;
        }
        result=0.0;
        for(i=0;i<max_ind;i++){
            result+=(double)myrank*(double)myrank;
        }
        // result = myrank;
        // MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&my_walkers[0], 1, MPI_WALKER, 0, 0, MPI_COMM_WORLD);
        work_counter++;
    }
    free(my_walkers);
    free(manager_walkers);
    fprintf(stderr, "myrank: %d, nwork: %d\n", myrank, work_counter);
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

    /*========================================================================*/
    /* make an MPI custom structure */
    int npars = (int)NPARS;
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
    if (myrank == 0) {
        manager(ntasks, MPI_WALKER);
    }
    else {
        worker(myrank, MPI_WALKER);
    }
    MPI_Type_free(&MPI_WALKER);
    MPI_Finalize();
}