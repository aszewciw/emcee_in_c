/*
    for docs see mca.h

    Copyright (C) 2012  Erin Sheldon

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "randn.h"
#include "mca.h"

struct mca_chain *mca_chain_new(size_t nwalkers,
                                size_t steps_per_walker,
                                size_t npars)
{
    struct mca_chain *self=calloc(1,sizeof(struct mca_chain));
    if (self==NULL) {
        fprintf(stderr,"Could not allocate struct mca_chain\n");
        exit(EXIT_FAILURE);
    }

    self->pars=calloc(nwalkers*steps_per_walker*npars,sizeof(double));
    if (self->pars==NULL) {
        fprintf(stderr,"Could not allocate mca_chain pars\n");
        exit(EXIT_FAILURE);
    }
    self->lnprob=calloc(nwalkers*steps_per_walker,sizeof(double));
    if (self->pars==NULL) {
        fprintf(stderr,"Could not allocate mca_chain lnprob\n");
        exit(EXIT_FAILURE);
    }
    self->accept=calloc(nwalkers*steps_per_walker,sizeof(int));
    if (self->pars==NULL) {
        fprintf(stderr,"Could not allocate mca_chain accept\n");
        exit(EXIT_FAILURE);
    }


    self->nwalkers=nwalkers;
    self->steps_per_walker=steps_per_walker;
    self->npars=npars;

    return self;

}
struct mca_chain *mca_chain_free(struct mca_chain *self)
{
    if (self) {
        free(self->pars);
        free(self->lnprob);
        free(self->accept);
        free(self);
    }
    return NULL;
}

int mca_chain_write_file(const struct mca_chain *self, const char *fname)
{
    FILE *fobj=fopen(fname,"w");
    if (fobj==NULL) {
        fprintf(stderr,"chain_writ error: Could not open file '%s'\n", fname);
        return 0;
    }

    mca_chain_write(self, fobj);
    fclose(fobj);

    return 1;
}
void mca_chain_write(const struct mca_chain *chain, FILE *stream)
{
    size_t nwalkers=MCA_CHAIN_NWALKERS(chain);
    size_t steps_per_walker=MCA_CHAIN_WNSTEPS(chain);
    size_t npars=MCA_CHAIN_NPARS(chain);

    size_t nsteps=nwalkers*steps_per_walker;

    fprintf(stream,"%lu %lu %lu\n", nwalkers, steps_per_walker, npars);
    for (size_t istep=0; istep<nsteps; istep++) {
        fprintf(stream,"%d %.16g ",
                MCA_CHAIN_ACCEPT(chain,istep),
                MCA_CHAIN_LNPROB(chain,istep));
        for (size_t ipar=0; ipar<npars; ipar++) {
            fprintf(stream,"%.16g ", MCA_CHAIN_PAR(chain, istep, ipar));
        }
        fprintf(stream,"\n");
    }
}

struct mca_chain *mca_chain_read(const char* fname)
{
    int status=1;
    struct mca_chain *chain=NULL;
    FILE *fptr=fopen(fname,"r");
    if (fptr==NULL) {
        fprintf(stderr,"chain_read error: Could not open file '%s'\n", fname);
        return NULL;
    }
    int nread=0;

    size_t nwalkers=0, steps_per_walker=0, npars=0;
    nread=fscanf(fptr,"%lu %lu %lu",
                 &nwalkers,&steps_per_walker,&npars);
    if (nread!=3) {
        fprintf(stderr,"Failed to read chain from %s",fname);
        status=0;
        goto _mca_chain_read_bail;
    }

    size_t nsteps_tot=nwalkers*steps_per_walker;
    chain=mca_chain_new(nwalkers,steps_per_walker,npars);

    for (size_t istep=0; istep<nsteps_tot; istep++) {
        nread=fscanf(
                fptr,
                "%d %lf",
                &MCA_CHAIN_ACCEPT(chain,istep),
                &MCA_CHAIN_LNPROB(chain,istep));
        if (nread != 2) {
            fprintf(stderr,"Failed to read chain from %s",fname);
            status=0;
            goto _mca_chain_read_bail;
        }

        nread=0;
        for (size_t ipar=0; ipar<npars; ipar++) {
            nread+=fscanf(fptr,"%lf", &MCA_CHAIN_PAR(chain, istep, ipar));
        }
        if (nread!=npars) {
            fprintf(stderr,"Failed to read chain from %s",fname);
            status=0;
            goto _mca_chain_read_bail;
        }
    }
_mca_chain_read_bail:
    if (status != 1) {
        // chain is null
        chain=mca_chain_free(chain);
    }
    fclose(fptr);

    return chain;
}


void mca_chain_plot(const struct mca_chain *self, const char *options)
{
    char cmd[256];
    char *name= tempnam(NULL,NULL);


    printf("writing temporary chain to: %s\n", name);

    FILE *fobj=fopen(name,"w");
    mca_chain_write(self, fobj);
    fclose(fobj);

    sprintf(cmd,"mca-plot %s %s", options, name);
    printf("%s\n",cmd);
    int ret=system(cmd);

    sprintf(cmd,"rm %s", name);
    printf("%s\n",cmd);
    ret=system(cmd);
    printf("ret: %d\n", ret);

    free(name);

}

struct mca_chain *mca_make_guess(double *centers,
                                 double *widths,
                                 size_t npars,
                                 size_t nwalkers)
{
    struct mca_chain *chain=mca_chain_new(nwalkers,1,npars);

    for (size_t ipar=0; ipar<npars; ipar++) {

        double center=centers[ipar];
        double width=widths[ipar];

        for (size_t iwalk=0; iwalk<nwalkers; iwalk++) {
            double val = center + width*(randu()-0.5)*2;
            MCA_CHAIN_WPAR(chain, iwalk, 0, ipar) = val;
        }
    }

    return chain;
}


struct mca_stats *mca_stats_new(size_t npars)
{
    struct mca_stats *self=calloc(1,sizeof(struct mca_stats));
    if (self==NULL) {
        fprintf(stderr,"Could not allocate struct mca_stats\n");
        exit(EXIT_FAILURE);
    }

    self->mean = calloc(npars, sizeof(double));
    if (self->mean == NULL) {
        fprintf(stderr,"Could not allocate mca_stats mean array\n");
        exit(EXIT_FAILURE);
    }
    self->cov = calloc(npars*npars, sizeof(double));
    if (self->cov == NULL) {
        fprintf(stderr,"Could not allocate mca_stats cov array\n");
        exit(EXIT_FAILURE);
    }

    self->npars=npars;
    return self;
}

struct mca_stats *mca_stats_free(struct mca_stats *self)
{
    if (self) {
        free(self->mean);
        free(self->cov);
        free(self);
    }
    return NULL;
}

struct mca_stats *mca_chain_stats(const struct mca_chain *chain)
{
    size_t npars = MCA_CHAIN_NPARS(chain);
    struct mca_stats *self=mca_stats_new(npars);

    if (!mca_chain_stats_fill(self,chain)) {
        self=mca_stats_free(self);
    }
    return self;
}

void mca_stats_clear(struct mca_stats *self)
{
    if (!self)
        return;

    memset(self->mean, 0, self->npars*sizeof(double));
    memset(self->cov, 0, self->npars*self->npars*sizeof(double));
    self->arate=0;

}
int mca_chain_stats_fill(
        struct mca_stats *self,
        const struct mca_chain *chain)
{
    size_t npars = MCA_CHAIN_NPARS(chain);
    size_t nsteps = MCA_CHAIN_NSTEPS(chain);
    double ival=0, jval=0;

    if (self->npars != npars) {
        fprintf(stderr,"mca_chain_stats error: Expected "
                       "npars %lu got %lu\n", npars, self->npars);
        return 0;
    }

    mca_stats_clear(self);

    double arate=0;
    for (size_t istep=0; istep<nsteps; istep++) {

        arate += MCA_CHAIN_ACCEPT(chain, istep);
        for (size_t ipar=0; ipar<npars; ipar++) {

            ival=MCA_CHAIN_PAR(chain,istep,ipar);
            self->mean[ipar] += ival;

            for (size_t jpar=ipar; jpar<npars; jpar++) {

                if (ipar==jpar) {
                    jval=ival;
                } else {
                    ival=MCA_CHAIN_PAR(chain,istep,jpar);
                }

                self->cov[ipar*npars + jpar] += ival*jval;
            }

        }
    }

    for (size_t ipar=0; ipar<npars; ipar++) {
        self->mean[ipar] /= nsteps;
    }

    for (size_t ipar=0; ipar<npars; ipar++) {
        double imean=self->mean[ipar];

        for (size_t jpar=ipar; jpar<npars; jpar++) {
            size_t index=ipar*npars + jpar;

            double jmean=self->mean[jpar];

            self->cov[index] /= nsteps;
            self->cov[index] -= imean*jmean;

            if (ipar!=jpar) {
                self->cov[jpar*npars + ipar] = self->cov[index];
            }
        }
    }

    self->arate = arate/nsteps;
    return 1;
}

void mca_stats_write_brief(struct mca_stats *self, FILE *stream)
{
    if (!self)
        return;
    size_t npars = MCA_STATS_NPARS(self);

    fprintf(stream,"%.16g\n", MCA_STATS_ARATE(self));
    for (size_t ipar=0; ipar< npars; ipar++) {
        double mn=MCA_STATS_MEAN(self,ipar);
        double var=MCA_STATS_COV(self,ipar,ipar);
        double err=sqrt(var);
        fprintf(stream,"%g +/- %g\n",mn,err);
    }
}

void mca_stats_write(struct mca_stats *self, FILE *stream)
{
    if (!self)
        return;
    size_t npars = MCA_STATS_NPARS(self);

    fprintf(stream,"%lu\n", npars);
    fprintf(stream,"%.16g\n", MCA_STATS_ARATE(self));
    for (size_t ipar=0; ipar< npars; ipar++) {
        double mn=MCA_STATS_MEAN(self,ipar);
        fprintf(stream,"%.16g ", mn);
    }
    fprintf(stream,"\n");
    for (size_t ipar=0; ipar< npars; ipar++) {
        for (size_t jpar=0; jpar< npars; jpar++) {
            double cov=MCA_STATS_COV(self,ipar,jpar);
            fprintf(stream,"%.16g ",cov);
        }
        fprintf(stream,"\n");
    }
}
// write space separated, no new line at all
void mca_stats_write_flat(struct mca_stats *self, FILE *stream)
{
    if (!self)
        return;
    size_t npars = MCA_STATS_NPARS(self);

    fprintf(stream,"%lu ", npars);
    fprintf(stream,"%.16g ", MCA_STATS_ARATE(self));
    for (size_t ipar=0; ipar< npars; ipar++) {
        double mn=MCA_STATS_MEAN(self,ipar);
        fprintf(stream,"%.16g ", mn);
    }
    for (size_t ipar=0; ipar< npars; ipar++) {
        for (size_t jpar=0; jpar< npars; jpar++) {
            double cov=MCA_STATS_COV(self,ipar,jpar);
            fprintf(stream,"%.16g ",cov);
        }
    }
}

static void copy_pars(const double *self, double *pars_dst, size_t npars)
{
    memcpy(pars_dst, self, npars*sizeof(double));
}


 /* copy the last step in the start chain to the first step
   in the chain */
static void set_start(struct mca_chain *chain,
                      const struct mca_chain *start,
                      double (*lnprob)(const double *, size_t, const void *),
                      const void *userdata)
{
    size_t steps_per=MCA_CHAIN_WNSTEPS(start);
    size_t nwalkers=MCA_CHAIN_NWALKERS(start);
    size_t npars=MCA_CHAIN_NPARS(start);

    for (size_t iwalker=0; iwalker<nwalkers; iwalker++) {
        MCA_CHAIN_WACCEPT(chain,iwalker,0) = 1;

        double *pars_start = MCA_CHAIN_WPARS(start, iwalker, (steps_per-1));
        double *pars = MCA_CHAIN_WPARS(chain, iwalker, 0);

        copy_pars(pars_start, pars, npars);

        MCA_CHAIN_WLNPROB(chain,iwalker,0) = (*lnprob)(pars,npars,userdata);
    }
}

static int mca_accept(double lnprob_old,
                      double lnprob_new,
                      int npars,
                      double z)
{
    double lnprob_diff = (npars - 1.)*log(z) + lnprob_new - lnprob_old;
    double r = randu();

    if (lnprob_diff > log(r)) {
        return 1;
    } else {
        return 0;
    }
}


static void mca_stretch_move(double z,
                             size_t ndim,
                             double *pars_old,
                             const double *pars_comp,
                             double *pars_new)
{
    for (size_t i=0; i<ndim; i++) {

        double val=pars_old[i];
        double cval=pars_comp[i];

        //newpars[i] = cval + z*(val-cval);
        pars_new[i] = cval - z*(cval-val);
    }
}

/*

   get a random long index in [0,n) from the *complement* of the input
   current value, i.e. such that index!=current

*/

unsigned int mca_rand_complement(unsigned int n)
{
    // unsigned int i=current;
    // while (i == current) {
    //     i = genrand_uint32_max(n);
    // }
    unsigned int i=genrand_uint32_max(n);
    return i;
}
/*
static void step_walker(struct mca_chain *chain,
                        double a,
                        double (*lnprob)(const double *, size_t, const void *),
                        const void *userdata,
                        size_t istep, size_t iwalker)
{
    size_t npars=MCA_CHAIN_NPARS(chain);
    size_t nwalkers=MCA_CHAIN_NWALKERS(chain);
    size_t prev_step=istep-1;

    const double *pars_old = MCA_CHAIN_WPARS(chain,iwalker,prev_step);
    double lnprob_old = MCA_CHAIN_WLNPROB(chain,iwalker,prev_step);

    double *pars_new=MCA_CHAIN_WPARS(chain,iwalker,istep);

    long cwalker=mca_rand_complement(iwalker,nwalkers);
    const double *cpars=MCA_CHAIN_WPARS(chain,cwalker,prev_step);

    // this over-writes our chain at this step
    double z = mca_rand_gofz(a);
    mca_stretch_move(a, z, pars_old, cpars, npars, pars_new);

    double lnprob_new = (*lnprob)(pars_new,npars,userdata);

    int accept = mca_accept(lnprob_old, lnprob_new, npars, z);
    MCA_CHAIN_WACCEPT(chain,iwalker,istep) = accept;

    if (accept) {
        // we already set pars_new above
        MCA_CHAIN_WLNPROB(chain,iwalker,istep) = lnprob_new;
    } else {
        // copy the older pars over what we put in above
        copy_pars(pars_old, pars_new, npars);
        MCA_CHAIN_WLNPROB(chain,iwalker,istep) = lnprob_old;
    }
}
*/
static void step_walker(struct mca_chain *chain,
                        struct mca_chain *comp_chain,
                        double a,
                        double (*lnprob)(const double *, size_t, const void *),
                        const void *userdata,
                        size_t istep, size_t iwalker)
{
    size_t npars=MCA_CHAIN_NPARS(chain);
    size_t nwalkers=MCA_CHAIN_NWALKERS(chain);
    double pars_old[npars];
    double pars_new[npars];
    for(int ipar=0;ipar<npars;ipar++){
        pars_old[ipar] = MCA_CHAIN_WPAR(chain,iwalker,0,ipar);
    }

    double lnprob_old = MCA_CHAIN_WLNPROB(chain,iwalker,0);

    long cwalker=mca_rand_complement(nwalkers);
    const double *cpars=MCA_CHAIN_WPARS(comp_chain,cwalker,0);

    double z = mca_rand_gofz(a);

    mca_stretch_move(z, npars, pars_old, cpars, pars_new);

    double lnprob_new = (*lnprob)(pars_new,npars,userdata);

    int accept = mca_accept(lnprob_old, lnprob_new, npars, z);
    MCA_CHAIN_WACCEPT(chain,iwalker,0) = accept;

    if (accept) {
        for(int ipar=0;ipar<npars;ipar++){
            MCA_CHAIN_WPAR(chain,iwalker,0,ipar) = pars_new[ipar];
        }
        MCA_CHAIN_WLNPROB(chain,iwalker,0) = lnprob_new;
    } else {
        MCA_CHAIN_WLNPROB(chain,iwalker,0) = lnprob_old;
    }
}

/*
void mca_run(struct mca_chain *chain,
             double a,
             const struct mca_chain *start,
             double (*lnprob)(const double *, size_t, const void *),
             const void *userdata)
{
    size_t nwalkers=MCA_CHAIN_NWALKERS(chain);
    size_t steps_per_walker=MCA_CHAIN_WNSTEPS(chain);

    set_start(chain, start, lnprob, userdata);

    for (size_t istep=1; istep<steps_per_walker; istep++) {
        for (size_t iwalker=0; iwalker<nwalkers; iwalker++) {
            step_walker(chain,a,lnprob,userdata,istep,iwalker);
        }
    }

}
*/

void chain_to_subchains(struct mca_chain *chain,
                        struct mca_chain *sub_chain1,
                        struct mca_chain *sub_chain2,
                        size_t istep)
{
    size_t npars=MCA_CHAIN_NPARS(chain);
    // size_t steps_per=MCA_CHAIN_WNSTEPS(chain);
    // size_t nwalkers=MCA_CHAIN_NWALKERS(chain);
    size_t nwalkers_over_two=MCA_CHAIN_NWALKERS(sub_chain1);

    for (size_t iwalker=0; iwalker<nwalkers_over_two; iwalker++) {

        size_t iwalker2=iwalker+nwalkers_over_two;

        MCA_CHAIN_WACCEPT(sub_chain1,iwalker,0) = MCA_CHAIN_WACCEPT(chain,iwalker,istep);
        MCA_CHAIN_WACCEPT(sub_chain2,iwalker,0) = MCA_CHAIN_WACCEPT(chain,iwalker2,istep);
        MCA_CHAIN_WLNPROB(sub_chain1,iwalker,0) = MCA_CHAIN_WLNPROB(chain,iwalker,istep);
        MCA_CHAIN_WLNPROB(sub_chain2,iwalker,0) = MCA_CHAIN_WLNPROB(chain,iwalker2,istep);
        for (size_t ipar=0; ipar<npars; ipar++){
            MCA_CHAIN_WPAR(sub_chain1,iwalker,0,ipar) = MCA_CHAIN_WPAR(chain,iwalker,istep,ipar);
            MCA_CHAIN_WPAR(sub_chain2,iwalker,0,ipar) = MCA_CHAIN_WPAR(chain,iwalker2,istep,ipar);
        }
    }
}

void subchains_to_chain(struct mca_chain *chain,
                        struct mca_chain *sub_chain1,
                        struct mca_chain *sub_chain2,
                        size_t istep)
{
    size_t npars=MCA_CHAIN_NPARS(chain);
    // size_t steps_per=MCA_CHAIN_WNSTEPS(chain);
    // size_t nwalkers=MCA_CHAIN_NWALKERS(chain);
    size_t nwalkers_over_two=MCA_CHAIN_NWALKERS(sub_chain1);

    for (size_t iwalker=0; iwalker<nwalkers_over_two; iwalker++) {

        size_t iwalker2=iwalker+nwalkers_over_two;

        MCA_CHAIN_WACCEPT(chain,iwalker,istep) = MCA_CHAIN_WACCEPT(sub_chain1,iwalker,0);
        MCA_CHAIN_WACCEPT(chain,iwalker2,istep) = MCA_CHAIN_WACCEPT(sub_chain2,iwalker,0);
        MCA_CHAIN_WLNPROB(chain,iwalker,istep) = MCA_CHAIN_WLNPROB(sub_chain1,iwalker,0);
        MCA_CHAIN_WLNPROB(chain,iwalker2,istep) = MCA_CHAIN_WLNPROB(sub_chain2,iwalker,0);
        for (size_t ipar=0; ipar<npars; ipar++){
            MCA_CHAIN_WPAR(chain,iwalker,istep,ipar) = MCA_CHAIN_WPAR(sub_chain1,iwalker,0,ipar);
            MCA_CHAIN_WPAR(chain,iwalker2,istep,ipar) = MCA_CHAIN_WPAR(sub_chain2,iwalker,0,ipar);
        }
    }
}

void mca_run(struct mca_chain *chain,
             double a,
             const struct mca_chain *start,
             double (*lnprob)(const double *, size_t, const void *),
             const void *userdata)
{
    size_t nwalkers=MCA_CHAIN_NWALKERS(chain);
    size_t steps_per_walker=MCA_CHAIN_WNSTEPS(chain);
    size_t npars=MCA_CHAIN_NPARS(chain);

    size_t nwalkers_over_two=nwalkers/2;
    struct mca_chain *sub_chain1=mca_chain_new(nwalkers_over_two,1,npars);
    struct mca_chain *sub_chain2=mca_chain_new(nwalkers_over_two,1,npars);

    fprintf(stderr, "Sub chains initialized\n");
    set_start(chain, start, lnprob, userdata);
    chain_to_subchains(chain, sub_chain1, sub_chain2, 0);

    fprintf(stderr, "Sub-chains copied to main chain\n");
    for (size_t istep=1; istep<steps_per_walker; istep++) {
        for (size_t iwalker=0; iwalker<nwalkers; iwalker++) {
            step_walker(sub_chain1,sub_chain2,a,lnprob,userdata,istep,iwalker);
            fprintf(stderr, "Stepped walker\n");
            step_walker(sub_chain2,sub_chain1,a,lnprob,userdata,istep,iwalker);
        }
        subchains_to_chain(chain, sub_chain1, sub_chain2, istep);
    }
    mca_chain_free(sub_chain1);
    mca_chain_free(sub_chain2);
}


double mca_rand_gofz(double a)
{
    // ( (a-1) rand + 1 )^2 / a;

    double z = (a - 1.)*randu() + 1.;

    z = z*z/a;

    return z;
}

/*
   Generate normal random numbers.

   Note we get two per run but I'm only using one.
*/

double mca_randn()
{
    double x1, x2, w, y1;//, y2;

    do {
        x1 = 2.*randu() - 1.0;
        x2 = 2.*randu() - 1.0;
        w = x1*x1 + x2*x2;
    } while ( w >= 1.0 );

    w = sqrt( (-2.*log( w ) ) / w );
    y1 = x1*w;
    //y2 = x2*w;
    return y1;
}

