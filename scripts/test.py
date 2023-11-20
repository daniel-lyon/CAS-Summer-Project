from astropy.io import fits



# file = '/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/results/NGC_247/GCs/kb231024_00050_icubes_cut.fits'
file = '/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/results/NGC_247/GCs/obj1/50.fits'

header = fits.getheader(file)
data = fits.getdata(file)

print(header)

y = data
x = header['WAVGOOD0'] + header['CD3_3']

import matplotlib.pyplot as plt

plt.plot(x,y)
plt.show()
