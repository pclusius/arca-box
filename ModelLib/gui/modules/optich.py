"""
=============================================================================
Copyright (C) 2021  Multi-Scale Modelling group
Institute for Atmospheric and Earth System Research (INAR), University of Helsinki
Contact information arca@helsinki.fi

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
=============================================================================
"""


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

def plot(addr):
    fixed = False
    c=['orangered','seagreen','deepskyblue','darkorchid']
    if os.path.exists(os.path.join(addr,'Changes.txt')):
        time, max_d_dpar,max_d_npar,max_d_vap,max_d_npd = np.genfromtxt(os.path.join(addr,'Changes.txt'), unpack=True)

        if len(time)<2:
            return 'The file was a stub.'
        f = open(os.path.join(addr,'Changes.txt'))
        x = f.readline()
        f.close()
        x = x.split()
        fig, axes = plt.subplots(2,1, figsize=(8,6))
        if np.sum([float(i) for i in x[5:]]) == 0: fixed = True
        for j,ax in enumerate(axes):
            if sum(abs(max_d_dpar))>0:ax.plot(time, max_d_dpar*100,c=c[0],lw=2,label='max_d_dpar')
            if sum(abs(max_d_npar))>0:ax.plot(time, max_d_npar*100,c=c[1],lw=2,label='max_d_npar')
            if sum(abs(max_d_vap) )>0:ax.plot(time, max_d_vap *100, c=c[2],lw=2,label='max_d_vap')
            if sum(abs(max_d_npd) )>0:ax.plot(time, max_d_npd *100, c=c[3],lw=2,label='max_d_ndep')
            if len(x)>9:
                for i in range(3):
                    if not fixed:
                        if i==2:
                            ax.axhline((-1)*float(x[5+2*i])*1e2, c=c[i],linestyle='--', label=x[i+1]+' min')
                            ax.axhline(float(x[5+2*i])*1e2, c=c[i],linestyle=':',label=x[i+1]+' max')
                        else:
                            # ax.axhline(float(x[4+2*i])*1e2, c=c[i],linestyle='--', label=x[i+1]+' min')
                            ax.axhline(float(x[5+2*i])*1e2, c=c[i],linestyle=':',label=x[i+1]+' max')

            if j==1:
                ax.axhline(0, c='k', alpha=0.7)
                ax.set_yscale('symlog')
            else:
                if fixed:
                    ax.set_title('Relative changes in ONE TIMESTEP (used fixed timestep)')
                else:
                    ax.set_title('Relative changes in ONE TIMESTEP (used varying timestep)')
            ax.set_xlabel('time (s)')
            ax.set_xlabel('time (s)')
            ax.set_ylabel('Largest relative change (%)')
            ax.legend()
            ax.grid()

        plt.tight_layout()
        fig.savefig(os.path.join(addr,'Changes.png'), dpi=150)

        plt.close()
        return 'Saved figure in the run directory'
    else:
        return 'Directory "'+addr+'" did not contain the necessary file "Changes.txt".'
