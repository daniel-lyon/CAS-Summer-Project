# make a plot of galaxy luminosity and redshift
lum = [0.1, 0.2, 0.3, 0.4, 0.5]
z = [0.01, 0.02, 0.03, 0.04, 0.05]

import matplotlib.pyplot as plt

plt.scatter(lum, z)
plt.xlabel('Luminosity')
plt.ylabel('Redshift')
plt.show()

# reverse the axes
plt.scatter(z, lum)
plt.xlabel('Redshift')
plt.ylabel('Luminosity')
plt.show()