#include "emcee.h"

typedef struct mydata {
    double *data;
    double *ivar;
} mydata;

double lnprob(const double *pars, int npars, const void *userdata)
{
    int idata,jdata,ndata;
    double chi2,diff_i,diff_j,lnprob,cinv_ij;
    const mydata *data = userdata;

    chi2=0;

    for (idata=0; idata<npars; idata++) {
        diff_i = data->data[idata] - pars[idata];
        for(jdata=0; jdata<npars; jdata++){
            diff_j = data->data[jdata] - pars[jdata];
            cinv_ij = data->ivar[idata*npars+jdata];
            chi2 += diff_i*diff_j*cinv_ij;
        }
    }

    lnprob = -0.5*chi2;

    return lnprob;
}

int main( int argc, char ** argv )
{
    int ipar, jpar, iwalker, npars, nwalkers, nsteps, resume;
    double a = 2.0;
    char data_fname[] = "data/means.dat";
    char icov_fname[] = "data/icov.dat";
    char guess_fname[] = "data/guesses.dat";
    FILE *file;
    mydata *gaussian_data;
    walker_pos *start_pos;

    MPI_Init(&argc, &argv);

    resume = 0;

    npars = 5;
    nwalkers = 250;
    nsteps = 4000;

    gaussian_data = calloc(1,sizeof(mydata));
    if (gaussian_data==NULL) {
        fprintf(stderr,"Could not allocate struct mydata\n");
        exit(EXIT_FAILURE);
    }
    gaussian_data->data = calloc(npars,sizeof(double));
    if (gaussian_data->data==NULL) {
        fprintf(stderr,"Could not allocate data within gaussian_data\n");
        exit(EXIT_FAILURE);
    }
    gaussian_data->ivar = calloc(npars*npars,sizeof(double));
    if (gaussian_data->ivar==NULL) {
        fprintf(stderr,"Could not allocate ivar within gaussian_data\n");
        exit(EXIT_FAILURE);
    }

    /* read in means */
    if((file=fopen(data_fname,"r"))==NULL){
        fprintf(stderr, "Error: Cannot open file %s\n", data_fname);
        exit(EXIT_FAILURE);
    }
    for( ipar=0; ipar<npars; ipar++ ){
        fscanf(file, "%lf", &gaussian_data->data[ipar]);
    }
    fclose(file);

    /* read in inverse covariance matrix */
    if((file=fopen(icov_fname,"r"))==NULL){
        fprintf(stderr, "Error: Cannot open file %s\n", icov_fname);
        exit(EXIT_FAILURE);
    }
    for( ipar=0; ipar<npars; ipar++ ){
        for(jpar=0; jpar<npars; jpar++){
            fscanf(file, "%lf", &gaussian_data->ivar[ipar*npars+jpar]);
        }
    }
    fclose(file);

    /* make space for initial position */
    start_pos = allocate_walkers(nwalkers,npars,0);
    /* read in guesses */
    if((file=fopen(guess_fname,"r"))==NULL){
        fprintf(stderr, "Error: Cannot open file %s\n", guess_fname);
        exit(EXIT_FAILURE);
    }
    for( iwalker=0; iwalker<nwalkers; iwalker++){
        for (ipar=0; ipar<npars; ipar++){
            fscanf(file, "%lf", &start_pos[iwalker].pars[ipar]);
        }
    }
    fclose(file);

    const char fname[] = "data/nd_gaussian_chain_cversion.dat";

    // fprintf(stderr, "Ready to start chain\n");
    // start_pos = make_guess(guess,ballsize,nwalkers,npars);
    run_chain(nwalkers, nsteps, npars, resume, a, start_pos,
              &lnprob, gaussian_data, fname);

    MPI_Barrier(MPI_COMM_WORLD);
    free_walkers(nwalkers,start_pos);
    free(gaussian_data->data);
    free(gaussian_data->ivar);
    free(gaussian_data);
    MPI_Finalize();

    return 0;
}