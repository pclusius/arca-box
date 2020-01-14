import numpy as np
import scipy as sc
import scipy.interpolate


'''
This script combines different files to create input for SMAC.
The first file in input vector defines the time column; other files are interpolated to that time.
'''
date = '180411'
data_folder = '../input'
filename    = date+'.dat'

root = '/home/pecl/malte/dMalte/Malte_in/Box/April2018/PC'+date

input = ([
'SMEAR_TEMP_168.dat',
'SMEAR_RH.dat',
'SMEAR_CS.dat',
'SMEAR_PRESS.dat',
'SA_ground.dat',
'NH3.dat'
])

counter = 0
header = ''

for file in input:
    if counter == 0:
        time,meas = np.genfromtxt(root+'/'+file, usecols=(0,1), unpack=True)
        for i,t in enumerate(time[1:]):
            if t>time[i]: exit

        output_array = np.zeros((i+1, len(input)+1))
        output_array[:,0] = time[0:i+1]
        output_array[:,1] = meas[0:i+1]
        counter = 2
        header = header + '%24s%28s'%('time', file)
    else:
        time_this,meas = np.genfromtxt(root+'/'+file, usecols=(0,1), unpack=True)
        meas_i = sc.interpolate.interp1d(time_this,meas,kind='linear', fill_value=0,bounds_error=True)
        output_array[:,counter] = meas_i(time[0:i+1])
        header = header + '%28s'%(file)
        counter +=1
np.savetxt(data_folder+'/'+filename, output_array, delimiter='  ', fmt="%1.20e", header=header)
# plt.plot(output_array[:-1,0],output_array[:-1,3])
# plt.show()
