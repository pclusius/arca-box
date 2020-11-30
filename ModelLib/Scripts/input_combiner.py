import numpy as np
import scipy as sc
import scipy.interpolate


'''
This script combines different files to create input for THEBOX.
The first file in input vector defines the time column; other files are interpolated to that time.
'''
dates = [
'2016-06-12',
'2017-04-08',
'2017-05-28',
'2017-07-15',
'2018-05-31',
'2016-10-02',
'2017-05-23',
'2017-06-03',
'2018-05-27',
'2019-07-31',
'2017-03-19',
'2017-03-28',
'2017-04-04',
'2018-03-30',
'2018-06-12',
]


# input = ([
# 'SMEAR_TEMP_168.dat',
# 'SMEAR_PRESS.dat',
# 'SMEAR_RH.dat',
# 'SMEAR_SWR.dat',
# 'SA_tower.dat',
# 'SA_ground.dat',
# 'NH3.dat',
# 'SMEAR_CS.dat',
# 'SMEAR_GAS_CO.dat',
# 'SMEAR_GAS_NO2.dat',
# 'SMEAR_GAS_NO.dat',
# 'SMEAR_GAS_O3.dat',
# 'SMEAR_GAS_SO2.dat',
# ])
input = ([
'acetaldehyde.dat',
'acetic_acid.dat',
'acetone.dat',
'benzene.dat',
'CO.dat',
'ethanol.dat',
'isoprene.dat',
'mek.dat',
'methanol.dat',
'mono.dat',
'monofrag.dat',
'NO.dat',
'NOX.dat',
'O3.dat',
'pres.dat',
'rh.dat',
'SO2.dat',
'SWR.dat',
'temp.dat',
'toluene.dat',
])

combofile = 'VOC_0042.dat'

multicolumn = False

if not multicolumn:
    for date in dates:
        data_folder = '/home/pecl/01-TUTKIMUS/Gradu/Casestudy/'+date+'/input'
        column = 0
        header = ''
        source_dir = data_folder#'/home/pecl/04-MALTE/dMalte/Malte_in/Box/April2018/PC'+date
        filename    = 'all'+date+'.dat'
        for file in input:
            if column == 0:
                time,meas = np.genfromtxt(source_dir+'/'+file, usecols=(0,1), unpack=True)
                for i,t in enumerate(time[1:]):
                    if t>time[i]: exit

                output_array = np.zeros((i+1, len(input)+1))
                output_array[:,0] = time[0:i+1]
                output_array[:,1] = meas[0:i+1]
                column = 2
                header = header + '%24s%28s'%('time', file)
            else:
                time_this,meas = np.genfromtxt(source_dir+'/'+file, usecols=(0,1), unpack=True)
                meas_i = sc.interpolate.interp1d(time_this,meas,kind='linear', fill_value=0,bounds_error=True)
                output_array[:,column] = meas_i(time[0:i+1])
                header = header + '%28s'%(file)
                column +=1
        np.savetxt(data_folder+'/'+filename, output_array, delimiter='  ', fmt="%1.20e", header=header)

if multicolumn:
    for date in dates:
        data_folder = '../../INOUT/PC_20'+date[0:2]+'-'+date[2:4]+'-'+date[4:]+'/input'
        column = 0
        header = ''
        source_dir = '/home/pecl/04-MALTE/dMalte/Malte_in/Box/April2018/PC'+date
        filename    = 'voc'+date+'.dat'
        data = np.genfromtxt(source_dir+'/'+combofile)
        nc = np.shape(data)[1]
        f = open(data_folder+'/'+filename, 'w')

        previous = 0.0
        for i, line in enumerate(data):
            if line[0]>previous:
                zzz = '%24.12e'*nc %tuple(line)
                f.write(zzz+'\n')
            previous = line[0]
        f.close()


# if multicolumn:
#     for date in dates:
#         data_folder = '../../INOUT/HYDE/PC_20'+date[0:2]+'-'+date[2:4]+'-'+date[4:]+'/input'
#         column = 0
#         header = ''
#         source_dir = '/home/pecl/04-MALTE/dMalte/Malte_in/Box/April2018/PC'+date
#         filename    = 'voc'+date+'.dat'
#         for file in input:
#             if column == 0:
#                 data = np.genfromtxt(source_dir+'/'+file, unpack=True)
#                 for i,t in enumerate(data[0:]):
#                     if t>data[0,i]: exit
#
#                 output_array = np.zeros((i+1, len(input)+1))
#                 output_array[:,0] = time[0:i+1]
#                 output_array[:,1] = meas[0:i+1]
#                 column = 2
#                 header = header + '%24s%28s'%('time', file)
#             else:
#                 time_this,meas = np.genfromtxt(source_dir+'/'+file, usecols=(0,1), unpack=True)
#                 meas_i = sc.interpolate.interp1d(time_this,meas,kind='linear', fill_value=0,bounds_error=True)
#                 output_array[:,column] = meas_i(time[0:i+1])
#                 header = header + '%28s'%(file)
#                 column +=1
#         np.savetxt(data_folder+'/'+filename, output_array, delimiter='  ', fmt="%1.20e", header=header)

# plt.plot(output_array[:-1,0],output_array[:-1,3])
# plt.show()
