/*
   test fitting a constant
*/
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "../mca.h"

struct mydata {
    size_t ndata;
    const double *data;
    double ivar; // same error for each
};

double lnprob(const double *pars, size_t npars, const void *userdata)
{
    double chi2=0, diff=0;
    const struct mydata *mydata = userdata;

    for (size_t i=0; i<mydata->ndata; i++) {
        diff = mydata->data[i]-pars[0];
        chi2 += diff*diff*mydata->ivar;
    }

    double lnprob = -0.5*chi2;

    return lnprob;
}

int main(int argc, char **argv)
{
    double a=2;
    const char fname[]="test/chain.dat";
    time_t tm;

    size_t ndata=100;
    size_t npars=1;
    size_t nwalkers=20;
    size_t burn_per_walker=200;
    size_t steps_per_walker=200;

    double truepars[1] = {1};
    double guess[1] = {0};
    double ballsize[1] = {0};
    double fracerr=0.1;

    double err=fracerr*truepars[0];

    fprintf(stderr,"nwalkers:  %lu\n", nwalkers);
    fprintf(stderr,"burn per:  %lu\n", burn_per_walker);
    fprintf(stderr,"steps per: %lu\n", steps_per_walker);
    fprintf(stderr,"npars:     %lu\n", npars);
    fprintf(stderr,"truth:     %.16g\n", truepars[0]);
    fprintf(stderr,"npoints:   %lu\n", ndata);
    fprintf(stderr,"err per:   %.16g\n", err);
    fprintf(stderr,"expect err on mean: %.16g\n", err/sqrt(ndata));


    // set up the random number generator
    (void) time(&tm);
    srand48((long) tm);

    // set up the data
    double *data = malloc(ndata*sizeof(double));
    for (size_t i=0; i<ndata; i++) {
        data[i] = truepars[0] + err*mca_randn();
    }

    guess[0] = truepars[0] + err*mca_randn();
    ballsize[0] = 1.0;

    struct mca_chain *guesses=mca_make_guess(guess, ballsize, npars, nwalkers);

    struct mydata mydata;
    mydata.ndata = ndata;
    mydata.data = (const double*) data;
    mydata.ivar = 1/(err*err);

    struct mca_chain *burnin_chain=mca_chain_new(nwalkers,burn_per_walker,npars);
    mca_run(burnin_chain, a, guesses, &lnprob, &mydata);

    struct mca_chain *chain=mca_chain_new(nwalkers,steps_per_walker,npars);
    mca_run(chain, a, burnin_chain, &lnprob, &mydata);

    struct mca_stats *stats=mca_chain_stats(chain);
    fprintf(stderr,"\nStats:\n");
    mca_stats_write_brief(stats,stderr);
    fprintf(stderr,"\nStats full:\n");
    mca_stats_write(stats,stderr);

    fprintf(stderr,"writing chain to %s\n", fname);
    mca_chain_write_file(chain, fname);

    fprintf(stderr,"reading and checking chain\n");
    struct mca_chain *tchain=mca_chain_read(fname);
    if (MCA_CHAIN_NSTEPS(chain)!=MCA_CHAIN_NSTEPS(tchain)) {
        fprintf(stderr,"expected steps: %lu, got %lu\n",
                MCA_CHAIN_NSTEPS(chain),MCA_CHAIN_NSTEPS(tchain));
        exit(EXIT_FAILURE);
    }
    for (size_t istep=0; istep<MCA_CHAIN_NSTEPS(chain); istep++) {
        if (MCA_CHAIN_ACCEPT(chain,istep) != MCA_CHAIN_ACCEPT(tchain,istep)) {
            fprintf(stderr,"step %lu expected accept: %d, got %d\n",
                    istep,
                    MCA_CHAIN_ACCEPT(chain,istep),
                    MCA_CHAIN_ACCEPT(tchain,istep));
            exit(EXIT_FAILURE);
        }

        for (size_t ipar=0; ipar<npars; ipar++) {
            double pardiff=
                fabs( MCA_CHAIN_PAR(chain,istep,ipar)-
                        MCA_CHAIN_PAR(tchain,istep,ipar) );
            pardiff /= MCA_CHAIN_PAR(chain,istep,ipar);
            if (pardiff > 1.e-15) {
                fprintf(stderr,"step %lu par %lu rel diff %lf\n",
                        istep,ipar,pardiff);
                exit(EXIT_FAILURE);
            }
        }

    }

    guesses=mca_chain_free(guesses);
    burnin_chain=mca_chain_free(burnin_chain);
    chain=mca_chain_free(chain);

    stats=mca_stats_free(stats);
    free(data);
}

