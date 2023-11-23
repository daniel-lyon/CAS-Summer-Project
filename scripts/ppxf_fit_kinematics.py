# Editted version of ppxf_example_kinematics_sauron.py written by michelle cappellari
# Aim is to fit template models to versions of the VCC1287 data
# Jonah Gannon
# October 2019

import glob
from os import path
from time import perf_counter as clock
import scipy.ndimage as ndi
from astropy.io import fits
from scipy import ndimage
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import ppxf as ppxf_package
from ppxf.ppxf import ppxf
import ppxf.ppxf_util as util
class HaltException(Exception): pass
################################################# USER INPUTS ##########################################################

fittable_file = ("/mnt/d/Ubuntu/230915/kb_t01_only/redux/kb230915_00019_icubed_cut_qspectrum.fits")

noise_file = ''

data_file_out = '/mnt/d/Ubuntu/virgo_megadwarfs/ppxf_fits/VCC1448_mosaic_mask11_median_241ppxf_fits.dat'

library_name = 'Coelho'

region = 'kcwi_new' # define a region to fit (default/all/no_mgb/blue/red/before_mgb/no_starting_hbeta)

z = 0.0044 # Initial redshift estimate of the galaxy
sigma_start = 33.

moments = 2 # gauss hermite moments to fit
degree  = 6 # additive legendre polynomial to fit
mdegree = 6 # multiplicative legendre polynomial to fit

cwave = 4500 # central wavelength of FWHM diff # makes almost no difference
FWHM_data = 2.53 # note this is changed to 2.52 for miles (1.02 or 0.524 depending on config)
# use 0.524 (measured) for bh3/medium
# use 1.02  (measured) for bh3/large
# use 2.52 for miles (same as templates)
# use 0.99 for BM/M
# 2.53 for Bl/medium (from website r val)
# use 5.055 for BL/Large (from website R val)
# use 5.72 for BL/L degraded

smoothing = 'internal' #internal to ppxf of external with ndi.gaussian
noise_source = 'constant' # file or constant source of errors

write_out = False
########################################################################################################################

##### code block from bron to import conroy models

def get_SSP_library_new_conroy_models(filepath):
    """
    Get the SSP library for the new Conroy models - these have a weird file structure, so they need their own function.
    Input:
        filepath: (string) the general filepath to all of the ssp files, eg. '/models/ssp_mist_c3k/SSP*.spec'
    Returns:
        ssp_lib
        ssp_data
        ssp_pop_variables
        ssp_lamrange
    """
    #load the library names
    ssp_lib = glob.glob(filepath)
    ssp_data = []
    ssp_ages = []
    ssp_metals = []
    #read in the first file
    for file in ssp_lib:
        with open(file, 'r') as fp:
            index = 0
            for count, line in enumerate(fp):
                #get rid of the comment lines
                if line[0] == "#":
                    continue
                #get the size of the arrays and create arrays to fill with population variables and spectra
                if count == 8:
                    data_size = np.array([int(s) for s in line[:-1].split()])
                    ssp_age = np.zeros((data_size[0]))
                    ssp_metal = np.full_like(ssp_age, float(file.split('Z')[1][:6]))
                    #spectra needs to have shape [nPixels, nAge] for pPXF
                    ssp_spectra = np.zeros((data_size[0], data_size[1]))
                #get the wavelength vector
                if count == 9:
                    lamdas = np.array([float(s) for s in line.split()])
                    lamdas_mask = (lamdas>3000)&(lamdas<7000)
                    lamdas = lamdas[lamdas_mask]
                #get all of the spectra and population variables, and fill the arrays
                if count > 9:
                    #make the line into an array
                    line_array = np.array([float(s) for s in line.split()])
                    #if it has four values, put it in the variables
                    if line_array.shape[0] == 4:
                        ssp_age[index] = (10**line_array[0])/(10**9)
                    #if it has the same length as expected for the spectra, put it in the spectra array
                    elif line_array.shape[0] == data_size[1]:
                        ssp_spectra[index, :] = line_array
                        index += 1
                    else:
                        print('Line {}: shape {} - line is not the right length!'.format(count, line_array.shape[0]))
            #use the wavelength mask on the spectra
            ssp_spectra = ssp_spectra[:, lamdas_mask]
            #add all of this to the total lists
            ssp_data.append(ssp_spectra)
            ssp_ages.append(ssp_age)
            ssp_metals.append(ssp_metal)
        fp.close()
    #find the wavelength range
    ssp_lamrange = np.array([lamdas[0], lamdas[-1]])
    #turn the data, ages, metals into arrays
    ssp_data = np.array(ssp_data)
    ssp_ages = np.array(ssp_ages)
    ssp_metals = np.array(ssp_metals)
    #reshape the spectra so that they have [nTemplates, nPixels]
    ssp_data = ssp_data.reshape((-1, ssp_data.shape[-1]))
    #reshape the population properties so they have [nPixels]
    ssp_ages = ssp_ages.reshape(-1)
    ssp_metals = ssp_metals.reshape(-1)
    return ssp_lib, ssp_data, ssp_ages, ssp_metals, ssp_lamrange
####



################################## Auto Populate based on input library name ###########################################
if library_name == 'Coelho':
    FWHM_tem = 0.254
    templates_loc = '/mnt/d/Ubuntu/Libraries/Coelho Libraries/2014/all/'

elif library_name == 'Coelho-4000':
    FWHM_tem = 0.254
    templates_loc = '/mnt/d/Ubuntu/Coelho Libraries/2014/4000/'

elif library_name == 'MILES':
    FWHM_tem = 2.51
    FWHM_gal = 2.52
    templates_loc = '/mnt/d/Ubuntu/Libraries/MILES/'

elif library_name == 'ELODIE':
    FWHM_tem = 0.508
    templates_loc = '/mnt/d/Ubuntu/ELODIE/Library/'

elif library_name == 'ELODIE-HiRes':
    FWHM_tem = 0.121
    templates_loc = '/media/gannonjs/Jonah Data/HiRes/'

elif library_name == 'ELODIE-Upsampled':
    FWHM_tem = 0.508
    templates_loc = '/mnt/d/Ubuntu/ELODIE/Library/upsampled/'

elif library_name == 'PEGASE':
    FWHM_tem = 0.508
    templates_loc = '/mnt/d/Ubuntu/PEGASE_HR/'

elif library_name == 'UVES':
    FWHM_tem = 0.0635
    templates_loc = '/media/gannonjs/Jonah Data/uves_pop/ktypes/'

elif library_name == 'M3':
    FWHM_tem = FWHM_data
    templates_loc = '/mnt/d/Ubuntu/M3/'

elif library_name == 'M3-Medium':
    FWHM_tem = FWHM_data
    templates_loc = '/mnt/d/Ubuntu/M3-Medium/'

elif library_name == 'M13':
    FWHM_tem = FWHM_data
    templates_loc = '/mnt/d/Ubuntu/M13/'

elif library_name == 'c3k':
    FWHM_tem = 0.508
    templates_loc = '/mnt/d/Ubuntu/c3k/SSP*.spec'
else:
    raise HaltException("Jonah Check your inputs you numpty")

########################################################################################################################

# Read a galaxy spectrum and define the wavelength range
#
hdu = fits.open(fittable_file)
gal_lin = hdu[0].data
h1 = hdu[0].header

wavelength = np.linspace(h1['CRVAL1'], h1['CRVAL1']+h1['CDELT1'] * (gal_lin.shape[0]-1), gal_lin.shape[0])

lamRange1 = h1['CRVAL1'] + np.array([0., h1['CDELT1']*(h1['NAXIS1'] - 1)])
#lamRange1 = [4823.75, 5312.25] # Undo for M3

galaxy, logLam1, velscale = util.log_rebin(lamRange1, gal_lin)
galaxy_normaliser = np.abs(np.median(galaxy))
galaxy = galaxy/galaxy_normaliser  # Normalize spectrum to avoid numerical issues

if noise_source == 'constant':
    noise = np.full_like(galaxy, 1)  # Assume constant noise per pixel here
elif noise_source == 'file':
    noise_lin = fits.open(noise_file)[0].data
    noise, logLam1, velscale = util.log_rebin(lamRange1, noise_lin)
    noise = noise/galaxy_normaliser
else:
    raise HaltException("Jonah define your noise")

# define a resolution

FWHM_gal = FWHM_data

# Uncomment below to see a spectrum after log rebinning
#plt.plot(np.exp(logLam1), galaxy)
#plt.plot(wavelength, galaxy)
#plt.legend(("Logarithmic", "Linear"))
#plt.vlines(4880, -50, 50)
#plt.show()

# plot the data with the fill
fig = plt.figure(54, figsize = (10, 6))
ax = plt.subplot(111)
ax.plot(np.exp(logLam1), galaxy)
ax.fill_between(np.exp(logLam1), galaxy - noise, galaxy+noise, color = 'lightgray', alpha = 0.75)
ax.set_xlabel("Wavelength [$\mathrm{\AA}$]")
ax.set_ylabel("Relative Flux")
#plt.show()

if library_name == 'c3k':
    ssp_file_list, ssp_data, ssp_ages, ssp_metals, lamRange2 = get_SSP_library_new_conroy_models(templates_loc)

    ssp = ssp_data[0]
    sspNew, logLam2, velscale_temp = util.log_rebin(lamRange2, ssp, velscale=velscale)
    templates = np.empty((sspNew.size, len(ssp_data)))

    FWHM_dif = np.sqrt(FWHM_gal ** 2 - FWHM_tem ** 2)

    if smoothing == 'external':
        raise HaltException("No external smoothing for C3K")

    for k in range(0,len(ssp_data)):
        ssp = ssp_data[k]
        sspNew, logLam2, velscale_temp = util.log_rebin(lamRange2, ssp, velscale=velscale)
        templates[:, k] = sspNew/np.median(sspNew)

else:
    vazdekis = sorted(glob.glob(templates_loc + '*.fits'))

    hdu = fits.open(vazdekis[0])

    if library_name == 'MILES':
        ssp = hdu[0].data[0]
    else:
        ssp = hdu[0].data
    h2 = hdu[0].header
    lamRange2 = h2['CRVAL1'] + np.array([0., h2['CDELT1']*(h2['NAXIS1'] - 1)])
    sspNew, logLam2, velscale_temp = util.log_rebin(lamRange2, ssp, velscale=velscale)
    templates = np.empty((sspNew.size, len(vazdekis)))

    FWHM_dif = np.sqrt(FWHM_gal**2 - FWHM_tem**2)
    sigma = FWHM_dif/2.355/h2['CDELT1']  # Sigma difference in pixels

    for j, file in enumerate(vazdekis):
        hdu = fits.open(file)
        if library_name == 'MILES':
            ssp = hdu[0].data[0]
        else:
            ssp = hdu[0].data
        #helps to remove NANs when problematic

        if library_name == 'ELODIE':
            whereNan = np.isnan(ssp)
            ssp[whereNan] = np.nanmedian(ssp)

        if library_name == 'ELODIE-Upsampled':
            whereNan = np.isnan(ssp)
            ssp[whereNan] = np.nanmedian(ssp)

        if library_name == 'ELODIE-HiRes':
            whereNan = np.isnan(ssp)
            ssp[whereNan] = np.nanmedian(ssp)

        # removes nans when there
        if smoothing == 'external':
            ssp = ndimage.gaussian_filter1d(ssp, sigma)

        sspNew, logLam2, velscale_temp = util.log_rebin(lamRange2, ssp, velscale=velscale)
        templates[:, j] = sspNew/np.median(sspNew)  # Normalizes templates

c = 299792.458
dv = (np.mean(logLam2[:1]) - logLam1[0])*c  # km/s

if region == 'all':
    goodPixels = np.arange(0, len(galaxy)-10, 1) # all
elif region == 'default':
    goodPixels = util.determine_goodpixels(logLam1, lamRange2, z)
elif region == 'no_mgb':
    goodPixels = np.hstack((np.arange(0, 1430, 1), np.arange(1550, len(galaxy), 1)))
elif region == 'blue':
    goodPixels = np.arange(0, round(len(galaxy)/2), 1)
elif region == 'kcwi_new':
    goodPixels = np.arange(200, 1850, 1)
elif region == 'lydia':
    goodPixels = np.arange(277, 2113, 1)
elif region == 'red':
    goodPixels = np.arange(round(len(galaxy)/2), len(galaxy)-25, 1)
elif region == 'before_mgb':
    goodPixels = np.arange(0, 1410, 1)
elif region == 'no_starting_hbeta':
    goodPixels = np.arange(175, len(galaxy), 1)
elif region == 'leo':
    goodPixels = np.arange(150, len(galaxy), 1)
elif region == 'lt5100':
    goodPixels = np.arange(0, 1128, 1)
elif region == 'lt5100-leo':
    goodPixels = np.arange(0, 1103, 1)
else:
    raise HaltException("Jonah you numpty define a fitting region!")


# Here the actual fit starts. The best fit is plotted on the screen.
# Gas emission lines are excluded from the pPXF fit using the GOODPIXELS keyword.
#
vel = c*np.log(1 + z)   # eq.(8) of Cappellari (2017)
t = clock()

sigma_diff  = c / (cwave/FWHM_dif) /2.354

print(fittable_file)

print(""
      "######################## FITTING ##############################"
      "")

if smoothing == 'internal':
    pp1 = ppxf(templates,
              galaxy,
              noise,
              velscale,
            goodpixels=goodPixels,
            plot=True,
            moments=moments, #leave this alone
            degree=degree,
            mdegree=mdegree,
            vsyst=dv,
            component=np.zeros_like(len(templates)),
            start=[vel, sigma_start],  # (km/s), starting guess for [V, sigma]
            sigma_diff = sigma_diff
            )

elif smoothing == 'external':
    pp1 = ppxf(templates,
              galaxy,
              noise,
              velscale,
            goodpixels=goodPixels,
            plot=True,
            moments=moments, #leave this alone
            degree=degree,
            mdegree=mdegree,
            vsyst=dv,
            start=[vel, sigma_start]  # (km/s), starting guess for [V, sigma]
            )
else:
    raise HaltException("Jonah define your smoothing")

print("Formal errors:")
print("     dV    dsigma   dh3      dh4")
print("".join("%8.2g" % f for f in pp1.error*np.sqrt(pp1.chi2)))

print('Elapsed time in pPXF: %.2f s' % (clock() - t))

plt.show()

############################ Write out some stuff ######################################################################
# note the wavelengths show below are just wrong because they are all in the weird log rebin units
if write_out == True:
    wavelength = np.exp(logLam1)
    galaxy = pp1.galaxy
    bestfit = pp1.bestfit

    writeable = np.column_stack((wavelength, galaxy, bestfit))
    np.savetxt(data_file_out, writeable, fmt="%s")

########################################################################################################################

#estimate SN from first fit
residuals = galaxy - pp1.bestfit
signal = np.mean((pp1.bestfit - pp1.apoly))#/pp1.mpoly)
noise = np.std(residuals)

SN_ratio = signal / noise

print("Finally the estimates spectra signal to noise is: %.2f" % SN_ratio)