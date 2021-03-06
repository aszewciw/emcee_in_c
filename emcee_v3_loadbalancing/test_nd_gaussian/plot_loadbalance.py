chain_fname = 'nd_gaussian_chain_cversion.dat'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    nprocs=4
    nwalkers=20
    nsteps=20

    walker_id = np.zeros((nsteps,nwalkers))
    for i in range(nsteps):
        walker_id[i,:] = np.arange(nwalkers)

    step,rank = np.genfromtxt(chain_fname,usecols=[0,1],unpack=True,comments='#')

    rank = rank.reshape((nsteps,nwalkers))

    plt.clf()
    plt.figure(1)
    plt_name = 'load_balance_test.png'
    cmap = matplotlib.colors.ListedColormap(['red', 'cyan', 'blue'])
    bounds = [1,2,3,4]
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    plt.scatter(step,walker_id,c=rank,cmap=cmap, s=100, norm=norm)
    plt.axhline(9.5,color='k')
    plt.axvline(0.5,color='k')
    plt.axhline(12.5,color='g',ls='--')
    plt.axhline(2.5,color='g',ls='--')
    plt.axis([-1,20,-1,20])
    plt.xlabel('step number')
    plt.ylabel('walker id')
    plt.title('Load balancing test, color-coded by process id')
    plt.savefig(plt_name)

if __name__ == '__main__':
    main()