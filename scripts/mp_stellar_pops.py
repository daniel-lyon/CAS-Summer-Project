# aim here is to create a parameter exploration for ppxf and stellar pops like we do for kinematics
# Jonah Gannon swinburne research fellow
# November 2022, works on ubuntu 22.04 lts

import os
from astropy.io import fits
import numpy as np
import ppxf.miles_util as miles_util
import ppxf as ppxf
import ppxf.ppxf_util as util
from multiprocessing import Pool
import warnings

warnings.filterwarnings(category=DeprecationWarning, action="ignore")
warnings.filterwarnings(category=RuntimeWarning, action="ignore")
warnings.filterwarnings(category=np.VisibleDeprecationWarning, action="ignore")

parent = os.getcwd()
model_dir = parent + '/MILES_BASTI_KU_baseFe/'
fittable_file = parent + '/results_GC/F8D1/obj1/mean_NCS.fits'
np_array_out = parent + '/F8D1_fit_inputs_no_restframe_all.npy'

region = "all" # define a region to fit (default/all/no_mgb/blue/red/before_mgb/lt5100)

#### FOR DF9 L_BL_4550 (out of focus)
FWHM_data = 4550 / 1800 # data resolution
start0=-108 ;  start1=60
mask_l=3586; mask_h=5575 # wavelength masks to apply at either end

n_balmer = 3  ;  n_forbidden = 3   # need to know

# loop over the polynomials
degree_min = 0 # additive legendre polynomial to fit
degree_max = 15
mdegree_min = 0 # multiplicative legendre polynomial to fit
mdegree_max = 15

########################################################################################################################

adegree = np.linspace(degree_min, degree_max, degree_max+1)
mdegree = np.linspace(mdegree_min, mdegree_max, mdegree_max+1)

###-----for all spectra the same
FWHM_tem = 2.5  # spectral resolution of MILES (change if using degraded models)
c = 299792.458  # Speed of light in kms-1
z = 0

### read if fits
hdu = fits.open(fittable_file)
t = hdu[0].data
th = hdu[0].header
##### Only use the wavelength range in common between galaxy and stellar library.
wave = th['CRVAL1'] + (np.arange(0., (th['NAXIS1'])) - 1) * th['CDELT1']
# wave = wave / (1 + zi)
lamRange = np.array([np.min(wave), np.max(wave)])
frac = wave[1] / wave[0]  # Constant lambda fraction per pixel
velscale = np.log(frac) * c  # Constant velocity scale in km/s per pixel

FWHM_dif = np.sqrt(FWHM_tem ** 2 - FWHM_data ** 2)
mask = (wave > mask_l) & (wave < mask_h)
flux = t[mask]
wave = wave[mask]
lamRange = np.array([np.min(wave), np.max(wave)])

galaxy, logwave, VelScale = util.log_rebin(lamRange, flux, velscale=None)
# print("Velscale = %s" % VelScale)
wave = np.exp(logwave)
# galaxy = galaxy/np.median(galaxy)            # Normalize spectrum to avoid numerical issues
noise = np.full_like(galaxy, 0.01)  # Assume constant noise per pixel here

#########################------------------- Setup templates -----------------------
pathname = model_dir + '/Mku1.30*baseFe*.fits'  # EMILES Basti_KU_baseFe
nmod = 516
miles = miles_util.miles(pathname, VelScale, FWHM_data, FWHM_tem)
reg_dim = miles.templates.shape[1:]
stars_templates = miles.templates.reshape(miles.templates.shape[0], -1)
lam_range_gal = np.array([np.min(wave), np.max(wave)]) / (1 + z)
gas_templates, gas_names, line_wave = \
util.emission_lines(miles.log_lam_temp, lam_range_gal, FWHM_tem)
templates = np.column_stack([stars_templates, gas_templates])
dv = c * (miles.log_lam_temp[0] - np.log(wave[0]))  # km/s
### .determine_goodpixels needs ln(gal_lambda) vector, minmax of template lin-lambda, z
lamRangeTemp = np.exp(np.array([np.min(miles.log_lam_temp), np.max(miles.log_lam_temp)]))

if region == 'all':
    goodPixels = np.arange(0, len(galaxy), 1) # all
elif region == "all_minus_10":
    goodPixels = np.arange(10, len(galaxy)-10, 1) # all
elif region == 'default':
    goodPixels = util.determine_goodpixels(logwave, lamRangeTemp, z)
elif region == 'no_mgb':
    goodPixels = np.hstack((np.arange(0, 1430, 1), np.arange(1550, len(galaxy), 1)))
elif region == 'blue':
    goodPixels = np.arange(0, round(len(galaxy)/2), 1)
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
    raise Exception("Jonah you numpty define a fitting region!")

storage_array = np.array([]) # first row to be the good fit, following rows to be the error fits with random masks

#########---------------- Fitting  -------------------------------------------
SigmaStart = 3 * VelScale  # if a good initial guess of the LOSVD parameters is not available
vel = c * np.log(1 + z)  # eq.(8) of Capp*n_teellari (2017)
start = [start0, start1]
n_temps = stars_templates.shape[1]
component = [0] * n_temps + [1] * (n_balmer + n_forbidden)  # will give us both balmer+forbidden together
gas_component = np.array(component) > 0  # gas_component=True for gas templates

moments = [2, 2]  # if balmer+forbiden
start = [start, start]
fixed = [[1, 1], [0, 0]]  # we will fix the Vr and sigma

def fit_and_save_results_parallel(params):
    i, deg_k, j, deg_p = params
    tracker = i * len(mdegree)+ j + 1
    print("\n"
            "######################## FIT %s of %s ##############################"
            "" % (str(tracker), str(len(adegree)*len(mdegree))))
    print("\n"
            "Degree: %s / mDegree: %s" % (
            str(i), str(j)))

    # Fit to get the kinematics and fix them
    pp = ppxf.ppxf.ppxf(templates, galaxy, noise, VelScale, start, goodpixels=goodPixels,
                        plot=False, moments=moments, degree=np.int(deg_k), mdegree=np.int(deg_p), vsyst=dv, lam=wave,
                        clean=False, regul=False, reg_dim=reg_dim, component=component,
                        gas_component=gas_component, gas_names=gas_names)  # ,velscale_ratio=velscale_ratio)

    noise1 = noise * np.sqrt(pp.chi2)

    # fit with fixed kinematics to get stellar pops
    pp1 = ppxf.ppxf.ppxf(templates, galaxy, noise1, VelScale, [pp.sol[0], pp.sol[1]], fixed=fixed, goodpixels=goodPixels,
                         plot=True, moments=moments, degree=np.int(deg_k), mdegree=np.int(deg_p), vsyst=dv, lam=wave,
                         clean=True, regul=1. / 0.1, reg_dim=reg_dim, component=component,
                         gas_component=gas_component, gas_names=gas_names)
    # plt.show()
    weights = pp1.weights[~gas_component]
    weights = weights.reshape(reg_dim) / weights.sum()
    mwt, mwz = miles.mean_age_metal(weights)
    snr = np.median(pp1.bestfit) / np.std(galaxy - pp1.bestfit)
    # print("SNR = %.2f" % (snr))
    collapsed_ages = np.sum(weights, axis=1)
    cumulative_growth = 1 - np.cumsum(collapsed_ages)
    result = np.array([i, j, mwz, mwt, snr, cumulative_growth])
    return result

if __name__ == "__main__":

    storage_array = np.array([])  # first row to be the good fit, following rows to be the error fits with random masks

    all_params = []
    for i, deg_k in enumerate(adegree):
        for j, deg_p in enumerate(mdegree):
            all_params.append([i, deg_k, j, deg_p])

    with Pool() as pool:
        # results = pool.map(fit_and_save_results_parallel, all_params)
        jobs = [pool.apply_async(fit_and_save_results_parallel, [param]) for param in all_params]
        results = [result.get() for result in jobs]

    storage_array = np.array(results)
    np.save(np_array_out, storage_array)
    print("All Done :)")