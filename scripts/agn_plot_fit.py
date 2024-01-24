import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# file_path1 = "/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/results_GC/NGC_247/5P/obj1/mean_NCS_smooth3.fits"
# file_path2 = "/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/results_GC/NGC_247/10P/obj1/mean_NCS_smooth3.fits"
# file_path3 = "/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/results_GC/IKN_GCC7/obj1/mean_NCS_smooth3.fits"
# names = ['5P', '10P', 'IKN']
# redshifts = [1.581, 0.909, 0.2335]

# mask_l = 3600
# mask_h = 5500

# files = [file_path1, file_path2, file_path3]
# n_files = len(files)
# fig, axes = plt.subplots(n_files, 1, figsize=(7, 5*n_files), sharex=True)

# for file, ax, name, z in zip(files, axes, names, redshifts):
#     hdu = fits.open(file)
#     flux = hdu[0].data
#     header = hdu[0].header
#     wave = header['CRVAL1'] + (np.arange(0., (header['NAXIS1'])) - 1) * header['CDELT1']
#     mask = (wave > mask_l) & (wave < mask_h)
#     wave = wave[mask]
#     flux = flux[mask]

#     # Plot the data
#     ax.plot(wave, flux, color="black", drawstyle="steps-mid")
#     ax.margins(x=0, y=0.3)
#     ax.tick_params(axis='both', which='major', labelsize=12)
#     ax.set_title(f'{name}, z = {z}', fontsize=12, y=0.875, va='top')

#     ax.axvline(x=1908.734 * (1 + z), color="grey", linestyle="dashed", linewidth=0.5) # CIII
#     ax.axvline(x=1857.4 * (1 + z), color="grey", linestyle="dashed", linewidth=0.5) # AlIII
#     ax.axvline(x=1549.48 * (1 + z), color="grey", linestyle="dashed", linewidth=0.5) # CIV
#     ax.axvline(x=1486 * (1 + z), color="grey", linestyle="dashed", linewidth=0.5) # NIV

#     ax.axvline(x=2799.117 * (1 + z), color="grey", linestyle="dashed", linewidth=0.5) # MgII

#     ax.axvline(x=4341.68 * (1 + z), color="grey", linestyle="dashed", linewidth=0.5) # Hg
#     ax.axvline(x=4102.89 * (1 + z), color="grey", linestyle="dashed", linewidth=0.5) # Hd
#     ax.axvline(x=3868.7 * (1 + z), color="grey", linestyle="dashed", linewidth=0.5) # NeIII
#     ax.axvline(x=3729.875 * (1 + z), color="grey", linestyle="dashed", linewidth=0.5) # OII
#     ax.axvline(x=3426.85 * (1 + z), color="grey", linestyle="dashed", linewidth=0.5) # NeVI

#     ax.text(x=(1908.734+10) * (1 + z), y=0.26, s="C III", fontsize=12, rotation='vertical', color="black")

#     if name != '10P':
#         ax.text(x=(1857.4+10) * (1 + z), y=0.26, s="Al III", fontsize=12, rotation='vertical', color="black")
#     ax.text(x=(1549.48+10) * (1 + z), y=0.26, s="C IV", fontsize=12, rotation='vertical', color="black")
#     ax.text(x=(1486+10) * (1 + z), y=0.26, s="N IV", fontsize=12, rotation='vertical', color="black")
    
#     if name != 'IKN':
#         ax.text(x=(2799.117+10) * (1 + z), y=0.26, s="Mg II", fontsize=12, rotation='vertical', color="black")

#     ax.text(x=(4341.68+10) * (1 + z), y=0.26, s="H\u03B3", fontsize=12, rotation='vertical', color="black")
#     ax.text(x=(4102.89+10) * (1 + z), y=0.26, s="H\u03B4", fontsize=12, rotation='vertical', color="black")
#     ax.text(x=(3868.7+10) * (1 + z), y=0.26, s="Ne III", fontsize=12, rotation='vertical', color="black")
#     ax.text(x=(3729.875+10) * (1 + z), y=0.26, s="O II", fontsize=12, rotation='vertical', color="black")
#     ax.text(x=(3426.85+10) * (1 + z), y=0.26, s="Ne VI", fontsize=12, rotation='vertical', color="black")

#     ax.set_xlim(wave[0], wave[-1])

# fig.text(0.02, 0.5, 'Relative Flux', va='center', rotation='vertical', fontsize=12)
# plt.xlabel("Wavelength $\AA$", fontsize=12)
# plt.subplots_adjust(hspace=0)
# plt.show()












file = "/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/results_GC/NGC_247/GCs2/obj1/spec.txt"
name = "NGC247 GC4"
v = 77

mask_l = 3600
mask_h = 5500

data = np.loadtxt(file)
x = data[:, 0]
galaxy = data[:, 1]
fit = data[:, 2]

# Plot the data
plt.figure(figsize=(10, 5))
plt.plot(x, galaxy, color="black", drawstyle="steps-mid")
plt.plot(x, fit, color="magenta", drawstyle="steps-mid")
plt.margins(x=0)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.title(f'{name}, v = {v} km/s', fontsize=12, y=0.95, va='top')
plt.xlabel("Wavelength $\AA$", fontsize=12)
plt.ylabel("Relative Flux", fontsize=12)
plt.show()
