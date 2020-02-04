import matplotlib.pyplot as plt

import numpy as np


plt.plot(np.linspace(0,1,100), np.sin(2*np.pi*np.linspace(0,1,100)))
plt.show()
