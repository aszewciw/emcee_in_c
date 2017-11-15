'''
code taken from http://dfm.io/emcee/current/user/quickstart/
'''
import emcee
import numpy as np

def lnprob(x, mu, icov):
    diff = x-mu
    return -np.dot(diff,np.dot(icov,diff))/2.0


chain_fname = 'nd_gaussian_chain_pyversion.dat'
ndim = 5

# ensure reproducibility for tests against python emcee
rseed = 10
np.random.seed(rseed)

means = np.random.rand(ndim)

cov = 0.5 - np.random.rand(ndim ** 2).reshape((ndim, ndim))
cov = np.triu(cov)
cov += cov.T - np.diag(cov.diagonal())
cov = np.dot(cov,cov)

icov = np.linalg.inv(cov)

nwalkers = 250
p0 = np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim))

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[means, icov])

lnprob(p, means, icov)

sampler.run_mcmc(p0, 4000)

np.savetxt(chain_fname,sampler.flatchain)