# a simple code that reads in fits files and outputs the barycentric correction
# jonah gannon phd student swinburne uni
# Jan 2020
# works on ubuntu 18.04 LTS


from barycentric import x_keckhelio
from astropy.io import fits
from astropy import wcs


############################################ USER INPUTS ###############################################################

# file_loc = '/mnt/d/Ubuntu/virgo_megadwarfs/VCC1448/yale_data/ocubes/pointing_b/kb210416_00085_ocubes.fits'
file_loc = '/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/data/Glob_Data_Cubes/NGC_243/GCs/kb231024_00053_icubes_cut.fits'

########################################################################################################################

header = fits.open(file_loc)[0].header

w = wcs.WCS(header)

ra_deg = w.wcs.crval[0]

dec_deg = w.wcs.crval[1]

Julian_Date = header["MJD"] + 2400000.5

bary_corr = x_keckhelio(ra = ra_deg, dec=dec_deg, jd = Julian_Date, obs = 'keck')

print("The Barycentric correction is: %.2f km/s" % bary_corr)
print("Be wary although - you must subtract this number")