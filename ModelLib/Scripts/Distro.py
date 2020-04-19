import numpy as np
import matplotlib.pyplot as plt

def g(sigma, mu):
  return np.exp(-( (d-mu)**2 / ( 2.0 * sigma**2 ) ) )

d = np.linspace(np.log10(1), np.log10(10000), 1000)

nucl = g(np.log10(3), np.log10(8))
aitk = g(np.log10(40), np.log10(100))
accu = g(np.log10(200), np.log10(400))

plt.semilogx(10**d, nucl)
plt.semilogx(10**d, aitk)
plt.semilogx(10**d, accu)
plt.semilogx(10**d, nucl + accu + aitk)
plt.show()
