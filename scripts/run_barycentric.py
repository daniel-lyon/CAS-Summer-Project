# a simple code that reads in fits files and outputs the barycentric correction
# jonah gannon phd student swinburne uni
# Jan 2020
# works on ubuntu 18.04 LTS


from barycentric import x_keckhelio
from astropy.io import fits
from astropy import wcs
import glob

########################################################################################################################

def barycentric_correction(file_directory, file_keyword):

    corrections = []
    file_list = sorted(glob.glob(file_directory + file_keyword))

    for file in file_list:

        header = fits.open(file)[0].header

        w = wcs.WCS(header)

        ra_deg = w.wcs.crval[0]

        dec_deg = w.wcs.crval[1]

        Julian_Date = header["MJD"] + 2400000.5

        bary_corr = x_keckhelio(ra = ra_deg, dec=dec_deg, jd = Julian_Date, obs = 'keck')[0]

        corrections.append(-bary_corr)

        print(f"The Barycentric correction is: {bary_corr} km/s")

    return corrections

if __name__ == '__main__':
    file_directory = '/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/results/NGC_247/GCs2/'
    bc = barycentric_correction(file_directory, file_keyword='*icubes*')
    print(bc)