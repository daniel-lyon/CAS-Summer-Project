# Created to apply barycentric corrections to data and stack them
# By Jonah Gannon
# July 2019


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import glob
import time
import spectres

start = time.time()

############################################# User Inputs ##############################################################

file_directory = '/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/outdata/NGC_247/'
file_keyword = '*_obj3.fits'

mean_out         = '/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/outdata/NGC_247/obj3_mean.fits'
mean_out_NCS     = '/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/outdata/NGC_247/obj3_mean_NCS.fits'
median_out       = '/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/outdata/NGC_247/obj3_median.fits'
stddev_out       = '/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/outdata/NGC_247/obj3_stddev.fits'
stddev_out_NCS   = '/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/outdata/NGC_247/obj3_stddev_NCS.fits'

weights = 'Equal' # apply a weighting if wanted to the stacking, make sure they are in the same order as things are read in

#be wary of the sign for below, by convention the barycentric correction is subtracted
bary_corr = np.array([0, 0, 0])
c = 299792 # km/s

write_out = True

########################################################################################################################

if weights == 'Equal' or 'equal':
    weights = np.ones_like(bary_corr)
else:
    pass

# Create File list
file_list = sorted(glob.glob(file_directory + file_keyword))
print(file_list)

spectra_data = []
spectra_header = []

# Import the spectra

for i, file in enumerate(file_list):
    spectra_data.append(fits.open(file)[0].data)
    spectra_header.append(fits.open(file)[0].header)

#Define the wavelength vector used based on the bounds of the imported UDG spectra

wavelength = np.linspace(spectra_header[0]['CRVAL1'], spectra_header[0]['CRVAL1']+spectra_header[0]['CDELT1'] * (spectra_data[0].__len__()-1), spectra_data[0].__len__())
print(wavelength)

bary_wavelength = wavelength[3:-4]

bary_corr_spectra_data = []
bary_corr_spectra_data_noContSub = []

for i, data in enumerate(spectra_data):
    bary_corr_spectra_dummy = spectres.spectres(bary_wavelength, wavelength - (wavelength * bary_corr[i]/c), spectra_data[i])
    bary_corr_spectra_dummy_median = bary_corr_spectra_dummy - np.median(bary_corr_spectra_dummy)
    bary_corr_spectra_data.append(bary_corr_spectra_dummy_median)
    bary_corr_spectra_data_noContSub.append(bary_corr_spectra_dummy)

mean_stacked_spectra = np.average(bary_corr_spectra_data, axis = 0, weights = weights)
mean_stacked_spectra_no_cont_sub = np.average(bary_corr_spectra_data_noContSub, axis = 0, weights = weights)

median_stacked_spectra = np.median(bary_corr_spectra_data, axis =0)

stddev_stacked_spectra = np.std(bary_corr_spectra_data, axis =0)

stddev_stacked_spectra_NCS = np.std(bary_corr_spectra_data_noContSub, axis =0)

plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'axes.linewidth': 3})
plt.rcParams.update({'xtick.major.width':3})
plt.rcParams.update({'ytick.major.width':3})

plt.rcParams.update({'xtick.minor.width':3})
plt.rcParams.update({'ytick.minor.width':3})

plt.rcParams.update({'xtick.minor.size':8})
plt.rcParams.update({'ytick.minor.size':8})

plt.rcParams.update({'xtick.major.size':10})
plt.rcParams.update({'ytick.major.size':10})

fig = plt.figure(11, figsize = (32, 16))
ax = plt.subplot(111)

for i, bary_data in enumerate(bary_corr_spectra_data):
    ax.plot(bary_wavelength, bary_data, lw = 2)

ax.set_xlabel('Wavelength [$\AA$]')
ax.set_ylabel("Relative Flux")

plt.show()

fig = plt.figure(1, figsize=(32, 16))
ax = plt.subplot(111)

ax.plot(bary_wavelength, median_stacked_spectra, lw=2)
ax.plot(bary_wavelength, mean_stacked_spectra, lw=2)

ax.fill_between(bary_wavelength, mean_stacked_spectra - stddev_stacked_spectra, mean_stacked_spectra + stddev_stacked_spectra, color = 'lightgray', alpha = 0.75)


ax.set_xlabel('Wavelength [$\AA$]')
ax.set_ylabel("Relative Flux")

#ax.set_xlim([wavelength[664], wavelength[704]])

plt.legend(("Median Stack", "Weighted Mean Stack", "1 $\sigma$ fill"))

plt.show()

######################################################################################################################
# write out the final mean spectra

if write_out == True:

    hdu = fits.PrimaryHDU(data=mean_stacked_spectra)
    hdu.header['CRVAL1'] = bary_wavelength[0]
    hdu.header['CDELT1'] = (bary_wavelength[bary_wavelength.size - 1] - bary_wavelength[0]) / (bary_wavelength.size - 1)
    hdu.header['CUNIT1'] = 'ANGSTROM'

    hdu.writeto(mean_out)
######################################################################################################################

######################################################################################################################
# write out the final mean spectra

if write_out == True:
    hdu = fits.PrimaryHDU(data=mean_stacked_spectra_no_cont_sub)
    hdu.header['CRVAL1'] = bary_wavelength[0]
    hdu.header['CDELT1'] = (bary_wavelength[bary_wavelength.size - 1] - bary_wavelength[0]) / (
                bary_wavelength.size - 1)
    hdu.header['CUNIT1'] = 'ANGSTROM'

    hdu.writeto(mean_out_NCS)
######################################################################################################################

######################################################################################################################
# write out the final median spectra

    hdu = fits.PrimaryHDU(data=median_stacked_spectra)
    hdu.header['CRVAL1'] = bary_wavelength[0]
    hdu.header['CDELT1'] = (bary_wavelength[bary_wavelength.size - 1] - bary_wavelength[0]) / (bary_wavelength.size - 1)
    hdu.header['CUNIT1'] = 'ANGSTROM'

    hdu.writeto(median_out)

# write out the final std deviation spectra

    hdu = fits.PrimaryHDU(data=stddev_stacked_spectra)
    hdu.header['CRVAL1'] = bary_wavelength[0]
    hdu.header['CDELT1'] = (bary_wavelength[bary_wavelength.size - 1] - bary_wavelength[0]) / (bary_wavelength.size - 1)
    hdu.header['CUNIT1'] = 'ANGSTROM'

    hdu.writeto(stddev_out)
    print("File Written Out")

# write out the final std deviation spectra

    hdu = fits.PrimaryHDU(data=stddev_stacked_spectra_NCS)
    hdu.header['CRVAL1'] = bary_wavelength[0]
    hdu.header['CDELT1'] = (bary_wavelength[bary_wavelength.size - 1] - bary_wavelength[0]) / (bary_wavelength.size - 1)
    hdu.header['CUNIT1'] = 'ANGSTROM'

    hdu.writeto(stddev_out_NCS)
    print("File Written Out")

else:
    print("Nothing Written Out")
######################################################################################################################

#Finish timing code
end = time.time()

runtime = end-start
print('Code Competed Successfully in, %.2f seconds' % runtime)

