'''
code taken from http://dfm.io/emcee/current/user/quickstart/
'''

import numpy as np

data_fname = 'means.dat'
icov_fname = 'icov.dat'
guesses_fname = 'guesses.dat'
ndim = 3

# ensure reproducibility for tests against python emcee
rseed = 10
np.random.seed(rseed)

means = np.random.rand(ndim)

cov = 0.5 - np.random.rand(ndim ** 2).reshape((ndim, ndim))
cov = np.triu(cov)
cov += cov.T - np.diag(cov.diagonal())
cov = np.dot(cov,cov)

icov = np.linalg.inv(cov)

nwalkers = 10
guesses = np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim))

np.savetxt(data_fname,means,fmt='%.12e')
np.savetxt(icov_fname,icov,fmt='%.12e')
np.savetxt(guesses_fname,guesses,fmt='%.12e')