#include <stdlib.h>
#include <stdio.h>
#include "randn.h"
#include "mca.h"


int main(int argc, char **argv)
{
    size_t ntrial=1000000;
    unsigned int nmax=100;
    double a=2;

    randn_seed();

    if (argc > 1) {
        ntrial=atol(argv[1]);
    }

    unsigned int maxval=0;

    printf("testing mca_rand_long with %ld trials\n", ntrial);
    // first test the basic random number generator
    for (size_t i=0; i<ntrial; i++) {
        unsigned int r = rand_uint_max(nmax);
        if (r < 0 || r >= nmax) {
            fprintf(stderr,
                    "Error, found rand %u outside of range: [0,%u)\n", r, nmax);
            exit(EXIT_FAILURE);
        }

        if (r > maxval) {
            maxval=r;
        }
    }
    printf("got maxval %u from range [0,%u)\n", maxval, nmax);

    // now as a complement
    printf("testing mca_rand_complement with %ld trials\n", ntrial);
    unsigned int current=3;
    for (size_t i=0; i<ntrial; i++) {
        unsigned int r = mca_rand_complement(current, nmax);
        if (r==current) {
            fprintf(stderr,
                    "Error, rand equal to the complement: [0,%u)\n", current);
            exit(EXIT_FAILURE);
        }
    }

    printf("all tests passed\n");

    printf("generating randoms from g(z), for a=2.  you'll have to test how\n"
           "that looks yourself.  See ./gofz-points.txt\n");
    FILE *fobj=fopen("gofz-points.txt","w");
    for (long i=0; i<ntrial; i++) {
        double z=mca_rand_gofz(a);
        fprintf(fobj,"%.16g\n", z);
    }
    fclose(fobj);
}

