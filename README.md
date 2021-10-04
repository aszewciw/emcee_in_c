# emcee_in_c
Parallel (MPI) C implementation of the affine-invariant Markov chain Monte Carlo ensemble sampler [`emcee`](https://github.com/dfm/emcee).

# Description
`emcee_in_c` is a parallel C implementation of the affine-invariant Markov chain Monte Carlo (MCMC)
ensemble sampler first proposed by [Goodman & Weare (2010)](https://cims.nyu.edu/~weare/papers/d13.pdf)
and implemented in Python by [Foreman-Mackey et al. (2013)](https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract).

# Installation
* To install this software, simply clone the repository using:
    git clone https://github.com/aszewciw/emcee_in_c.git
* Alternatively, `emcee_in_c` can be set up as a submodule within another repository. The process for doing this is general and is described well in [this link](https://git-scm.com/book/en/v2/Git-Tools-Submodules).
* Use of this software requires an installation of MPI.

# How to Use
Forthcoming...

# Notes to the user
This software is free to use under the terms of the MIT License.
`emcee_in_c` does not have a dedicated release publication.
We ask that anyone making use of `emcee_in_c` please cite Szewciw et al. 2021 (in prep), where we introduce `emcee_in_c`.

We developed this C-implementation due to software limitations concerning the use of `python-mpi` on the Texas Advanced Computing Center's Stampede2 supercomputer.
As such, unless the user faces a similar issue, we do not recommend use of `emcee_in_c` over `emcee`, which has far superior development, documentation, and utilities.
Although `emcee_in_c` was primarily developed for personal use, I will attempt to respond to any issues raised by users of this software.

# Acknowledgements
We thank [Erin Sheldon](https://github.com/esheldon) for providing a non-MPI C implementation of `emcee` which we used as a starting point to develop `emcee_in_c`.