try:
    import numpy as np
except:
    print('Failed to load NumPy')
try:
    import matplotlib.pyplot as plt
except:
    print('Failed to load matplotlib')

import sys
import os

if len(sys.argv)>1:
    if os.path.exists(sys.argv[1]):
        os.chdir(sys.argv[1])

if os.path.exists('optimChanges.txt'):
    time, max_d_vap, max_d_npar, max_d_dpar = np.genfromtxt('optimChanges.txt', unpack=True)
    plt.plot(time, max_d_vap*100, label='max_d_vap')
    plt.plot(time, max_d_npar*100, label='max_d_npar')
    plt.plot(time, max_d_dpar*100, label='max_d_dpar')
    plt.axhline(0.5, c='k')
    plt.axhline(3, c='k')
    plt.axhline(-3,linestyle=':', c='k')
    plt.gca().set_yscale('symlog')
    # plt.ylim(-100,10)

    plt.legend()
    plt.show()
else:
    print('Directory did not contain the correct file')
