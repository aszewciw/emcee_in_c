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

    // fprintf(stderr, "data 0 is %lf\n", mydata->data[0]);
    // fprintf(stderr, "data 1 is %lf\n", mydata->data[1]);
    for (i=0; i<mydata->ndata; i++) {
        diff = mydata->data[i]-pars[0];
        chi2 += diff*diff*mydata->ivar;
    }

    lnprob = -0.5*chi2;

    return lnprob;
}

int main( int argc, char ** argv )
{
    int ndata=100;
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