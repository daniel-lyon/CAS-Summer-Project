# Used to crop the kcwi data cubes to the useful wavelength region
# By Jonah Gannon
# Nov 2019
# Swinburne PHD Student
# Tested on used on python 3/Lunix 18.04

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from math import ceil
from math import floor
import glob
import sys
import os.path as path

#Time the Code
import time
start = time.time()

############################################# USER INPUTS ##############################################################
cube_file_loc = '/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/data/Glob_Data_Cubes/NGC_243/GCs/'
cube_file_key = '*kb*'

append_text = '_cut' #what you append to the outfile

x_low = 7 # lowest and highest xvalue of the spaxel to be INCLUDED in the rectangle
x_high = 26

y_low = 17 # lowest and highest yvalue of the spaxel to be INCLUDED in the rectangle
y_high = 81

# rectangle allowed will then be xlow, ylow to xhigh, yhigh along the diagonal

write_out = True

########################################################################################################################

files_to_cut = sorted(glob.glob(cube_file_loc+cube_file_key))
print(files_to_cut)

# Want to calculate the best wavelength range to crop all of the cubes
wl_low = []
wl_high = []
begin_wl = []

for i, cube_file in enumerate(files_to_cut):
    cube = fits.open(cube_file)
    cube_data = cube[0].data
    pixels_per_AA = 1 / cube[0].header['CD3_3']

    wl_low.append(int((ceil(cube[0].header['WAVGOOD0']*pixels_per_AA)/pixels_per_AA - cube[0].header['CRVAL3']) / cube[0].header['CD3_3'])+2)
    wl_high.append(int((floor(cube[0].header['WAVGOOD1']*pixels_per_AA)/pixels_per_AA - cube[0].header['CRVAL3']) / cube[0].header['CD3_3']-2))

    begin_wl.append(cube[0].header['CRVAL3'])

    print(cube_file, begin_wl[i]) # untag for a useful diagnostic thingo

    if begin_wl[0] != begin_wl[i]:
        print("PROBLEMO - Wavelengths don't start in the same for the cubes!")
        sys.exit()

# So now take the best of those
wl_highest_low = np.max(wl_low)
wl_lowest_high = np.min(wl_high)
print("The beginning and ending wavelength cuts for the cube are %.2f and %.2f respectively" % (wl_highest_low, wl_lowest_high))

#some plotting stuffs
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
###

for i, cube_file in enumerate(files_to_cut):
    # Import Cube
    cube = fits.open(cube_file)

    cube_data = cube[0].data
    cube_data_cut = cube_data[wl_highest_low:wl_lowest_high, y_low-1:y_high, x_low-1:x_high]
    # import fits header
    cube_header = cube[0].header
    cube_header['CRVAL3'] = begin_wl[i] + wl_highest_low * cube[0].header['CD3_3']
    # cube_header["CRPIX3"] -= wl_highest_low - 2 # removed as it seems better to just update the starting pixel
    cube_header["CRPIX1"] = 1 + (x_high-x_low)/2
    cube_header["CRPIX2"] = 3 + (y_high-y_low)/2

    # Make a before/after plot
    fig = plt.figure(i, figsize=(16, 18))
    ax = plt.subplot(121)
    ax.imshow(np.median(cube_data, axis = 0), origin = 'lower')
    ax.set_title('Before Cropping')
    ax = plt.subplot(122)
    ax.imshow(np.median(cube_data_cut, axis = 0), origin = 'lower')
    ax.set_title('After Cropping')
    plt.suptitle('File is %s' % cube_file, y=1)
    plt.show()
    time.sleep(1)

    # now create a file name
    file_name = path.splitext(cube_file)
    file_name_final = file_name[0] + append_text + file_name[1]

    ######################################################################################################################
    # write out the final cube
    if write_out == True:
        hdu = fits.PrimaryHDU(data=cube_data_cut)
        hdu.header = cube_header
        hdu.writeto(file_name_final, overwrite=False)
        print("Files Written out")
    else:
        print("Shes done boss but nothing was written out")
    ######################################################################################################################
#Finish timing code
end = time.time()

runtime = end-start
print('Code Competed Successfully in, %.2f seconds' % runtime)