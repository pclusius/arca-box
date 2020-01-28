import numpy as np
import scipy as sc
import scipy.interpolate


'''
This script combines different files to create input for THEBOX.
The first file in input vector defines the time column; other files are interpolated to that time.
'''
dates = [
'180401',
'180402',
'180403',
'180404',
'180405',
'180405',
'180407',
'180408',
'180409',
'180410',
'180411',
'180412',
'180413',
'180414',
'180415',
'180416',
'180417',
'180418',
'180419',
'180420',
'180421',
'180422',
'180423',
'180424',
'180425',
'180426',
'180427',
'180428',
'180429',
'180430',
]



input = ([
'SMEAR_TEMP_168.dat',
'SMEAR_PRESS.dat',
'SMEAR_RH.dat',
'SMEAR_SWR.dat',
'SA_tower.dat',
'SA_ground.dat',
'NH3.dat',
'SMEAR_CS.dat',
'SMEAR_GAS_CO.dat',
'SMEAR_GAS_NO2.dat',
'SMEAR_GAS_NO.dat',
'SMEAR_GAS_O3.dat',
'SMEAR_GAS_SO2.dat',
])
# 'VOC_0042.dat',

for date in dates:
    data_folder = '../../INOUT/HYDE/PC_20'+date[0:2]+'-'+date[2:4]+'-'+date[4:]+'/input'
    counter = 0
    header = ''
    source_dir = '/home/pecl/04-MALTE/dMalte/Malte_in/Box/April2018/PC'+date
    filename    = 'env'+date+'.dat'
    for file in input:
        if counter == 0:
            time,meas = np.genfromtxt(source_dir+'/'+file, usecols=(0,1), unpack=True)
            for i,t in enumerate(time[1:]):
                if t>time[i]: exit

            output_array = np.zeros((i+1, len(input)+1))
            output_array[:,0] = time[0:i+1]
            output_array[:,1] = meas[0:i+1]
            counter = 2
            header = header + '%24s%28s'%('time', file)
        else:
            time_this,meas = np.genfromtxt(source_dir+'/'+file, usecols=(0,1), unpack=True)
            meas_i = sc.interpolate.interp1d(time_this,meas,kind='linear', fill_value=0,bounds_error=True)
            output_array[:,counter] = meas_i(time[0:i+1])
            header = header + '%28s'%(file)
            counter +=1
    np.savetxt(data_folder+'/'+filename, output_array, delimiter='  ', fmt="%1.20e", header=header)
# plt.plot(output_array[:-1,0],output_array[:-1,3])
# plt.show()
