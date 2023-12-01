# !/usr/bin/env python
##############################################################################
import sys
from os.path import expanduser
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import time
from astropy.io import ascii
import ppxf as ppxf
import ppxf.ppxf_util as util
import ppxf.miles_util as miles_util
import progressbar
class HaltException(Exception): pass

#################################################### User Inputs #######################################################

fittable_file = '/mnt/d/Ubuntu/DF9_Letter/ocubes_ngc1052_df9/kb220129_00061_ocubes_cut_galaxy_no_nucleus_qspectrumSS.fits'  ## M_BH_5080
model_dir = '/mnt/d/Ubuntu/Libraries/MILES_SSP/MILES_BASTI_KU_baseFe/'

data_file_out = '/home/gannonjs/PycharmProjects/Research_Codes/DF9_Letter/data/kb220129_00061_ocubes_cut_galaxy_no_nucleus_qspectrumSS_ppxf_fit_288.dat'

write_out = True

region = "all" # define a region to fit (default/all/no_mgb/blue/red/before_mgb/lt5100)

#### FOR DF9 L_BL_4550 (out of focus)
zi = 0.0056 #redshift guess
FWHM_data =  5.72 # data resolution
start0=0 ;  start1=10
mask_l=3554; mask_h=5574 # wavelength masks to apply at either end
deg_k=8 ; deg_p=8; reg_val=0.1 # final polynomial / multiplicitive / regularisations to apply
n_balmer = 3  ;  n_forbidden = 3   # need to know

########################################################################################################################

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
wave = wave / (1 + zi)
lamRange = np.array([np.min(wave), np.max(wave)])
frac = wave[1] / wave[0]  # Constant lambda fraction per pixel
velscale = np.log(frac) * c  # Constant velocity scale in km/s per pixel

FWHM_dif = np.sqrt(FWHM_tem ** 2 - FWHM_data ** 2)
mask = (wave > mask_l) & (wave < mask_h)
flux = t[mask]
wave = wave[mask]
lamRange = np.array([np.min(wave), np.max(wave)])

galaxy, logwave, VelScale = util.log_rebin(lamRange, flux, velscale=None)
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
    raise HaltException("Jonah you numpty define a fitting region!")


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

####----------------  FOR KINEM DEG   -------------------------------------------
# deg_kinem = np.asarray([0,1,2,3,4,5,6,7,8,9,10,15,20,25])  #series of degrees to check
# for i in range(len(deg_kinem)):
#    if i ==0:
#        vr   = np.zeros(len(deg_kinem))
#        evr  = np.zeros(len(deg_kinem))
#        sig  = np.zeros(len(deg_kinem))
#        esig = np.zeros(len(deg_kinem))
#
#    pp = ppxf.ppxf.ppxf(templates, galaxy, noise, VelScale, start, goodpixels=goodPixels,
#          plot=True, moments=moments, degree=deg_kinem[i], mdegree=-1, vsyst=dv, lam=wave,
#          clean=False, regul=False, reg_dim=reg_dim, component=component,
#          gas_component=gas_component, gas_names=gas_names)#,velscale_ratio=velscale_ratio)
#
#    noise1 = noise*np.sqrt(pp.chi2)
#    vr[i]  = pp.sol[0][0]  ; evr[i]  = pp.error[0][0]*np.sqrt(pp.chi2)
#    sig[i] = pp.sol[0][1]  ; esig[i] = pp.error[0][1]*np.sqrt(pp.chi2)
#
#    plt.clf()
#    pp.plot()
#    plt.show()
#    time.sleep(5)
#    plt.close()
#    # plt.savefig('find_deg_reg/ppxffit_finddegkin_R84_BLsh2_deg'+str(format(i, '04'))+'.png')
#    data=[deg_kinem, vr, evr, sig, esig]
#    # ascii.write(data, 'find_deg_reg/ppxftable_r0_degkin_R84_BLsh2.txt', overwrite=True)
#    print('The deg_kin was %.1f'% deg_kinem[i])
#
# fig = plt.figure(1, figsize=(12, 8))
# ax1 = plt.subplot(211)
# ax2 = plt.subplot(212)
#
# ax1.errorbar(data[0], data[1], yerr = data[2])
# ax1.set_ylabel("Vr")
#
# ax2.errorbar(data[0], data[3], yerr = data[4])
# ax2.set_ylabel("$\sigma$")
# ax2.set_xlabel("Deg")
# plt.show()

########### Select the kinematic degree above then comment out above block and un-uncomment block below


########### Once we select the degree for the kinematics
pp = ppxf.ppxf.ppxf(templates, galaxy, noise, VelScale, start, goodpixels=goodPixels,
      plot=False, moments=moments, degree=deg_k, mdegree=-1, vsyst=dv, lam=wave,
      clean=False, regul=False, reg_dim=reg_dim, component=component,
      gas_component=gas_component, gas_names=gas_names)  # ,velscale_ratio=velscale_ratio)

noise1 = noise * np.sqrt(pp.chi2)

# then we need to select the multiplicative degree for the kinematics

#####----------------------------------------------------------------------
## ---------------  FOR POP DEG   -----------------------------------------
# --------------------------------------------------------------------------
# deg_pop = np.asarray([1,10])  #series of degrees to check
# deg_pop = np.asarray([1,2,3,4,5,6,7,8,9,10,15,20,25])  #series of degrees to check
# for i in range(len(deg_pop)):
#    if i ==0:
#        mwt  = np.zeros(len(deg_pop))
#        mwz  = np.zeros(len(deg_pop))
#        lwt  = np.zeros(len(deg_pop))
#        lwz  = np.zeros(len(deg_pop))
#        snr  = np.zeros(len(deg_pop))
#    pp1 = ppxf.ppxf.ppxf(templates, galaxy, noise1, VelScale, [pp.sol[0], pp.sol[1]], fixed=fixed, goodpixels=goodPixels,
#          plot=False, moments=moments, degree=deg_k, mdegree=deg_pop[i], vsyst=dv, lam=wave,
#          clean=True, regul=False, reg_dim=reg_dim, component=component,
#          gas_component=gas_component, gas_names=gas_names)#,velscale_ratio=velscale_ratio)
#
#    weights = pp1.weights[~gas_component]
#    weights = weights.reshape(reg_dim)/weights.sum()
#    mwt[i], mwz[i] = miles.mean_age_metal(weights)
#    snr[i] = np.median(pp1.bestfit) / np.std(galaxy - pp1.bestfit)
#
#    plt.clf()
#    plt.subplot(211)
#    pp1.plot()
#    plt.subplot(212)
#    miles.plot(weights)
#    plt.tight_layout()
#    plt.show()
#    time.sleep(5)
#    plt.close()
#    # plt.savefig('find_deg_reg/ppxffit_finddegpop_R84_BLsh1_deg'+str(format(i, '04'))+'.png')
#    data=[deg_pop, mwt, mwz, snr]
#    # ascii.write(data, 'find_deg_reg/ppxftable_r0_degpop_R84_BLsh1.txt', overwrite=True)
#    print('The deg_pop was %.1f'% deg_pop[i])
#
# fig = plt.figure(2, figsize=(12, 8))
# ax1 = plt.subplot(311)
# ax2 = plt.subplot(312)
# ax3 = plt.subplot(313)
#
# ax1.errorbar(data[0], data[1])
# ax1.set_ylabel("Mw Age")
#
# ax2.errorbar(data[0], data[2])
# ax2.set_ylabel("Mw [Z/H]")
#
# ax3.errorbar(data[0], data[3])
# ax3.set_ylabel("SNR")
#
# ax2.set_xlabel("Deg")
# plt.show()

##### Select the thing from above than uncomment below and start playing around with regularisations

##### and we fix the pop deg to look for the regularization
pp1 = ppxf.ppxf.ppxf(templates, galaxy, noise1, VelScale, [pp.sol[0], pp.sol[1]], fixed=fixed, goodpixels=goodPixels,
       plot=True, moments=moments, degree=deg_k, mdegree=deg_p, vsyst=dv, lam=wave,
       clean=True, regul=False, reg_dim=reg_dim, component=component,
       gas_component=gas_component, gas_names=gas_names)  # ,velscale_ratio=velscale_ratio)
plt.clf()
plt.subplot(211)
pp1.plot()
plt.subplot(212)
weights = pp1.weights[~gas_component]
weights = weights.reshape(reg_dim)/weights.sum()
miles.plot(weights)
plt.tight_layout()
plt.suptitle("Pops Degree Results")
plt.show()


mwt, mwz = miles.mean_age_metal(weights)
snr = np.median(pp1.bestfit) / np.std(galaxy - pp1.bestfit)
plt.show()

########----------------
######### B) and now we change the regul parameter until reduced and absolute are the same
#########------------------------------------------------------------------------------
# regul_err = np.asarray([1.0, 0.1, 0.01, 0.001, 0.0001, 0.00001])
# # regul_err = np.asarray([0.000009, 0.000008,0.000007,0.000006,0.000005,0.000004,0.000003,0.000002])
# # regul_err = np.asarray([0.00009, 0.00008,0.00007,0.00006,0.00005,0.00004,0.00003,0.00002])
# # regul_err = np.asarray([0.009, 0.008,0.007,0.006,0.005,0.004,0.003,0.002])
# # regul_err = np.asarray([0.09, 0.08,0.07,0.06,0.05,0.04,0.03,0.02])
# # regul_err = np.asarray([0.9, 0.8,0.7,0.6,0.5,0.4,0.3,0.2])
#
# for i in range(len(regul_err)):
#     if i == 0:
#         mwt = np.zeros(len(regul_err))
#         mwz = np.zeros(len(regul_err))
#         lwt = np.zeros(len(regul_err))
#         lwz = np.zeros(len(regul_err))
#         snr = np.zeros(len(regul_err))
#         red_chi = np.zeros(len(regul_err))
#         reg_chi = np.zeros(len(regul_err))
#         des_chi = np.zeros(len(regul_err))
#         cur_chi = np.zeros(len(regul_err))
#
#     pp2 = ppxf.ppxf.ppxf(templates, galaxy, noise1, VelScale, [pp.sol[0], pp.sol[1]], fixed=fixed, goodpixels=goodPixels,
#                plot=False,
#                moments=moments, degree=deg_k, mdegree=deg_p, vsyst=dv, lam=wave, clean=True,
#                regul=1. / regul_err[i], reg_dim=reg_dim,
#                component=component, gas_component=gas_component,
#                gas_names=gas_names)  # ,velscale_ratio=velscale_ratio)
#
#     weights = pp2.weights[~gas_component]
#     weights = weights.reshape(reg_dim) / weights.sum()
#     mwt[i], mwz[i] = miles.mean_age_metal(weights)
#     snr[i] = np.median(pp2.bestfit) / np.std(galaxy - pp2.bestfit)
#     red_chi[i] = pp2.chi2
#     reg_chi[i] = pp2.chi2 * goodPixels.size
#     des_chi[i] = np.sqrt(2 * goodPixels.size)
#     cur_chi[i] = (pp2.chi2 - pp1.chi2) * goodPixels.size
#     print('The desired chi was %.8g' % des_chi[i])
#     print('The current chi was %.8g' % cur_chi[i])
#
#     plt.clf()
#     plt.subplot(211)
#     pp2.plot()
#     plt.subplot(212)
#     miles.plot(weights)
#     plt.tight_layout()
#     plt.suptitle("The Regularisation_err was %.5f" % regul_err[i])
#     plt.show()
#     time.sleep(5)
#     plt.close()
#
#     data = [regul_err, red_chi, reg_chi, des_chi, cur_chi, mwt, mwz, snr]
#
#     print('The regul_error was %.8f' % regul_err[i])


############# Finally do the fit with everything you just worked out
pp2 = ppxf.ppxf.ppxf(templates, galaxy, noise1, VelScale, [pp.sol[0],pp.sol[1]], fixed=fixed, goodpixels=goodPixels,plot=False,
   moments=moments, degree=deg_k, mdegree=deg_p, vsyst=dv, lam=wave, clean=True, regul=1./reg_val, reg_dim=reg_dim,
   component=component,gas_component=gas_component, gas_names=gas_names)#,velscale_ratio=velscale_ratio)

weights = pp2.weights[~gas_component]
weights = weights.reshape(reg_dim)/weights.sum()
mwt,mwz = miles.mean_age_metal(weights)
snr = np.median(pp2.bestfit) / np.std(galaxy - pp2.bestfit)
plt.clf()
plt.subplot(211)
pp2.plot()
plt.subplot(212)
miles.plot(weights)
plt.tight_layout()
plt.suptitle("Final Results")
plt.show()


#
######------------------------------------------------------------------------------

############################ Write out some stuff ######################################################################
# note the wavelengths show below are just wrong because they are all in the weird log rebin units
if write_out == True:
    wavelength = np.exp(logwave)
    galaxy = pp1.galaxy
    bestfit = pp1.bestfit

    writeable = np.column_stack((wavelength, galaxy, bestfit))
    np.savetxt(data_file_out, writeable, fmt="%s")

########################################################################################################################

print("Its Done")