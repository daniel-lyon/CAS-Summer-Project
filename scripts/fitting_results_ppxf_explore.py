# aim is to plot the results of the fitting of the spectra
# Jonah Gannon Swinburne Postdoc
# November 2022
# Works on ubuntu 22.04lts

import os
import numpy as np
import matplotlib.pyplot as plt
import corner

""" Very good """
parent = os.getcwd()
# nucleus_loc = parent + '/test/GC1_all.npy'
# nucleus_loc = parent + '/test/GC1_default_lit_inputs.npy'
# nucleus_loc = parent + '/test/DDO190.npy'
# nucleus_loc = parent + '/test/M31_B336.npy'
# nucleus_loc = parent + '/test/F8D1.npy'
# nucleus_loc = parent + '/results_SP/NGC247 GCs GC1.npy'
# nucleus_loc = '/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/results_SP/NGC247 GCs GC2.npy'
# nucleus_loc = '/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/results_SP/NGC247 GCs GC3.npy'
# nucleus_loc = '/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/results_SP/NGC247 GCs2 GC1.npy'
# nucleus_loc = '/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/results_SP/NGC247 GCs2 GC2.npy'
# nucleus_loc = '/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/results_SP/NGC247 GCs2 GC3.npy'
# nucleus_loc = '/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/results_SP/M31_H12.npy'
# nucleus_loc = '/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/results_SP/M31_PANDAS_41.npy'
# nucleus_loc = '/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/results_SP/Sextans_A_GC1.npy'

""" Questionable """
# nucleus_loc = parent + '/test/M31_B336_fit_inputs.npy'
# nucleus_loc = parent + '/test/M31_B336_fit_inputs_no_restframe.npy'
# nucleus_loc = parent + '/test/M31_B336_LIT_inputs.npy'
# nucleus_loc = parent + '/results_SP/DDO190.npy'
# nucleus_loc = parent + '/results_SP/F8D1.npy'
# nucleus_loc = parent + '/results_SP/M31_B336.npy'

# nucleus_loc = parent + '/test/GC1_fit_inputs_no_restframe_default.npy'
# nucleus_loc = parent + '/test/GC2_fit_inputs_no_restframe_all.npy'
# nucleus_loc = parent + '/test/GC3_fit_inputs_no_restframe_default.npy'
# nucleus_loc = parent + '/test/GC4_fit_inputs_no_restframe_all.npy'
# nucleus_loc = parent + '/test/GC5_fit_inputs_no_restframe_all.npy'
nucleus_loc = parent + '/test/GC6_fit_inputs_no_restframe_all.npy'

exclude_nans = False # When true excludes all failed nan values from the data

miles_ages = np.array([0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0])
nucleus_data = np.load(nucleus_loc, allow_pickle=True)
# nucleus_data = np.reshape(nucleus_data, (np.int(len(nucleus_data)/6), 6))
nucleus_data[:, 3] = 10** nucleus_data[:, 3] / 10**9

if exclude_nans == True:
    full_nucleus_data = []
    for i in nucleus_data:
        if np.isnan(i[2]) == False:
            full_nucleus_data.append(i)
        else:
            print("nan excluded")
    nucleus_data = np.array(full_nucleus_data)

cornerfig = corner.corner(nucleus_data[:, 2:5], quantiles=[0.16, 0.5, 0.84], labels = ["[Z/H] [dex]]", "Age [Gyr]", "S/N"],show_titles=True)
plt.suptitle("For the Nucleus")
plt.show()

print("Note Jonah the S/N is being effected by the final 3 pixels of the fitting, removing them drastically increases S/N")