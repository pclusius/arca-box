#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 13:24:51 2019

@author: carltonx
"""

import os
import sys
import numpy as np
import pandas as pd
import pylab as pl
import matplotlib.pyplot as plt

path = '/home/carltonx/Box-model/PC180330/'

input_temp =np.loadtxt(path+'SMEAR_TEMP_168.dat')

print(input_temp.shape)


plt.plot(input_temp)
plt.show()
