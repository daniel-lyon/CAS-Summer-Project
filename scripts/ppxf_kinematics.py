#!/usr/bin/env python
##############################################################################
#
# Usage example for the procedure PPXF, which implements the
# Penalized Pixel-Fitting (pPXF) method originally described in
# Cappellari M., & Emsellem E., 2004, PASP, 116, 138
#     http://adsabs.harvard.edu/abs/2004PASP..116..138C
# and upgraded in Cappellari M., 2017, MNRAS, 466, 798
#     http://adsabs.harvard.edu/abs/2017MNRAS.466..798C
#
# The example also shows how to include a library of templates
# and how to mask gas emission lines if present.

##############################################################################

from os import path
from time import perf_counter as clock

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# from ppxf.ppxf import ppxf
# import ppxf.ppxf_util as util
# import ppxf.sps_util as lib

from ppxf.ppxf import ppxf
import ppxf.ppxf_util as util
import ppxf.miles_util as lib

def ppxf_kinematics(file, fwhm_gal, degree=4, wav_min=0, wav_max=-1, fit='default', moments=2, mdegree=6):

    # Read a galaxy spectrum and define the wavelength range
    hdu = fits.open(file)
    gal_lin = hdu[0].data
    h1 = hdu[0].header
    lamRange1 = h1['CRVAL1'] + np.array([0., h1['CDELT1']*(h1['NAXIS1'] - 1)])
    
    # Use these lines if your spectrum is at low-z (z<0.01
    redshift_0 = 0                  # Ignore cosmological redshift for local galaxies
    redshift = 0.000515               # Initial redshift estimate of the galaxy

    galaxy, ln_lam1, velscale = util.log_rebin(lamRange1, gal_lin)
    galaxy = galaxy/np.median(galaxy)  # Normalize spectrum to avoid numerical issues
    noise = np.full_like(galaxy, 0.0047)           # Assume constant noise per pixel here

    # Setup stellar templates
    sps_name = 'emiles'

    # The templates span a much larger wavelength range. To save some
    # computation time, I truncate the spectra to a similar but slightly larger
    # range than the galaxy.
    lam_range_temp = [lamRange1[0]/1.02, lamRange1[1]*1.02]

    # Read SPS models file from my GitHub if not already in the ppxf package dir.
    # The SPS model files are also available here https://github.com/micappe/ppxf_data
    basename = f"spectra_{sps_name}_9.0.npz"
    ppxf_dir = path.dirname(path.realpath(util.__file__))
    filename = path.join(ppxf_dir, 'sps_models', basename)

    sps = lib.sps_lib(filename, velscale, fwhm_gal, wave_range=lam_range_temp)

    # Compute a mask for gas emission lines
    if fit == 'default':
        goodPixels = util.determine_goodpixels(ln_lam1, lam_range_temp, redshift)[wav_min:wav_max]
    elif fit == 'all':
        goodPixels = np.arange(0, len(galaxy)-10, 1)[wav_min:wav_max] # all

    # Here the actual fit starts. The best fit is plotted on the screen. Gas
    # emission lines are excluded from the pPXF fit using the GOODPIXELS
    # keyword.
    c = 299792.458
    vel = c*np.log(1 + redshift)   # eq.(8) of Cappellari (2017, MNRAS)
    start = [vel, 200.]  # (km/s), starting guess for [V, sigma]
    t = clock()

    pp = ppxf(sps.templates, galaxy, noise, velscale, start,
              goodpixels=goodPixels, plot=True, 
              moments=moments, degree=degree, mdegree=mdegree,
              lam=np.exp(ln_lam1),
              lam_temp=sps.lam_temp, component=np.zeros_like(len(sps.templates)))
    
    residuals = galaxy[:len(goodPixels)] - pp.bestfit[:len(goodPixels)]
    signal = np.median((pp.bestfit[:len(goodPixels)] - pp.apoly[:len(goodPixels)]))#/pp1.mpoly)
    noise = np.std(residuals[residuals<np.percentile(residuals, 90)])

    SN_ratio = signal / noise

    # The updated best-fitting redshift is given by the following
    # lines (using equations 5 of Cappellari 2022, arXiv, C22)
    errors = pp.error*np.sqrt(pp.chi2)  # Assume the fit is good chi2/DOF=1
    redshift_fit = (1 + redshift_0)*np.exp(pp.sol[0]/c) - 1  # eq. (5c) C22
    redshift_err = (1 + redshift_fit)*errors[0]/c            # eq. (5d) C22

    print("Formal errors:")
    print("     dV    dsigma   dh3      dh4")
    print("".join("%8.2g" % f for f in errors))
    print('Elapsed time in pPXF: %.2f s' % (clock() - t))
    prec = int(1 - np.floor(np.log10(redshift_err)))  # two digits of uncertainty
    print(f"Best-fitting redshift z = {redshift_fit:#.{prec}f} "
          f"+/- {redshift_err:#.{prec}f}")
    print(f"Signal-to-noise ratio per pixel: {SN_ratio:.2f}")
    print(signal, noise)
    plt.show()
    return pp.sol[0], pp.sol[1], redshift_fit, redshift_err


##############################################################################

if __name__ == '__main__':

    file = '/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/results/Globs/Sextans_A_GC1/obj1/mean_NCS.fits'
    fwhm_gal = 4925 / 4000
    wavcut = -1
    degree = 6 # legendre polynomial degree
    ppxf_kinematics(file, fwhm_gal, degree, wavcut, moments=2, mdegree=6, fit='all')