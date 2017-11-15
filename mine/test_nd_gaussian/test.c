#include "emcee.h"

struct mydata {
    double data[NPARS];
    double ivar[NPARS][NPARS]; // same error for each
};

double lnprob(const double *pars, int npars, const void *userdata)
{
    int idata,jdata,ndata;
    double chi2,diff_i,diff_j,lnprob,cinv_ij;
    const struct mydata *mydata = userdata;

    ndata=(int)NPARS;

    chi2=0;

    for (idata=0; idata<ndata; idata++) {
        diff_i = mydata->data[idata] - pars[idata];
        for(jdata=0; jdata<ndata; jdata++){
            diff_j = mydata->data[jdata] - pars[jdata];
            cinv_ij = mydata->ivar[idata][jdata];
            chi2 += diff_i*diff_j*cinv_ij;
        }
    }

    lnprob = -0.5*chi2;

    return lnprob;
}

int main( int argc, char ** argv )
{
    double a=2.0;
    char data_fname[]="means.dat";
    char icov_fname[]="icov.dat";
    char guess_fname[]="guesses.dat";
    FILE *file;
    struct mydata *gaussian_data;
    int npars, nwalkers, ipar, jpar, iwalker;
    walker_pos *start_pos;

    npars = (int)NPARS;
    nwalkers = (int)NWALKERS;

    gaussian_data = calloc(1,sizeof(mydata));
    if (gaussian_data==NULL) {
        fprintf(stderr,"Could not allocate struct mydata\n");
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
            fscanf(file, "%lf", &gaussian_data->ivar[ipar][jpar]);
        }
    }
    fclose(file);

    /* make space for initial position */
    start_pos = allocate_walkers(nwalkers);
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

    const char fname[]="nd_gaussian_chain_cversion.dat";

    // start_pos = make_guess(guess,ballsize,nwalkers,npars);
    run_chain(&argc, &argv, start_pos, a, &lnprob, gaussian_data, fname);
    free(gaussian_data);
    free(start_pos);

    return 0;
}