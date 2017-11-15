'''
code taken from http://dfm.io/emcee/current/user/quickstart/
'''

data_fname = 'means.dat'
icov_fname = 'icov.dat'
guesses_fname = 'guesses.dat'
chain_fname = 'nd_gaussian_chain_cversion.dat'
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import corner

def main():

    steps = np.genfromtxt(chain_fname,usecols=[2,3,4,5,6],skip_header=1,comments='#')
    truths = np.genfromtxt(data_fname)

    signif_levels = np.array([1.0,2.0,3.0])
    levels = 1.0 - np.exp(-0.5*signif_levels**2)
    quantiles = [0.16,0.84]

    plt.clf()
    plt.figure(1)
    plt_name = 'nd_gaussian_chain_cversion_contours.png'
    fig = corner.corner(steps,levels=levels,color='r',quantiles=quantiles,
                        plot_density=False,plot_datapoints=False,truths=truths,
                        fill_contours=True)
    plt.savefig(plt_name)



if __name__ == '__main__':
    main()