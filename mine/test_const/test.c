#include "emcee.h"

typedef struct mydata {
    int ndata;
    const double *data;
    double ivar; // same error for each
} mydata;

double lnprob(const double *pars, int npars, const void *userdata)
{
    int i;
    double chi2,diff,lnprob;
    const mydata *data = userdata;

    chi2=0;
    diff=0;

    // fprintf(stderr, "data 0 is %lf\n", mydata->data[0]);
    // fprintf(stderr, "data 1 is %lf\n", mydata->data[1]);
    for (i=0; i<data->ndata; i++) {
        diff = data->data[i]-pars[0];
        chi2 += diff*diff*data->ivar;
    }

    lnprob = -0.5*chi2;

    return lnprob;
}

int main( int argc, char ** argv )
{
    int ndata=10;
    int i, npars, nwalkers;
    double a=2.0;

    double truepars[1] = {1};
    double guess[1] = {0};
    double ballsize[1] = {0};
    double fracerr=0.1;
    walker_pos *start_pos;

    double err=fracerr*truepars[0];

    double *data_val = malloc(ndata*sizeof(double));
    for(i=0; i<ndata; i++) {
        data_val[i] = truepars[0] + err*mca_randn();
    }

    nwalkers = (int)NWALKERS;
    npars = (int)NPARS;

    mydata data;
    data.ndata = ndata;
    data.data = (const double*) data_val;
    data.ivar = 1/(err*err);

    guess[0] = truepars[0] + err*mca_randn();
    ballsize[0] = 1.0;

    start_pos = make_guess(guess, ballsize, nwalkers, npars);

    const char fname[]="chain.dat";
    run_chain(&argc, &argv, start_pos, a, &lnprob, &data, fname);
    free(start_pos);
    free(data_val);

    return 0;
}