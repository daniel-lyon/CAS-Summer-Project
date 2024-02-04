import os
from matplotlib import patches
import numpy as np
import matplotlib.pyplot as plt
parent = os.getcwd()



""" NGC 247"""
# file_path1 = parent + "/results_GC/NGC_247/GCs/obj1/spec.txt"
# file_path2 = parent + "/results_GC/NGC_247/GCs/obj2/spec.txt"
# file_path3 = parent + "/results_GC/NGC_247/GCs/obj3/spec.txt"
# file_path4 = parent + "/results_GC/NGC_247/GCs2/obj1/spec.txt"
# file_path5 = parent + "/results_GC/NGC_247/GCs2/obj2/spec.txt"
# file_path6 = parent + "/results_GC/NGC_247/GCs2/obj3/spec.txt"

# files = [file_path1, file_path2, file_path3, file_path4, file_path5, file_path6]
# names = ["GC1", "GC2", "GC3", "GC4", "GC5", "GC6"]
# redshifts = [0.000278, 0.0001811, 0.000137, 0.000256, 0.000286, 0.000305]


""" Other candidates """
file_path1 = parent + "/results_GC/DDO190/obj1/spec.txt"
file_path2 = parent + "/results_GC/F8D1/obj1/spec.txt"
file_path3 = parent + "/results_GC/M31_B336/obj1/spec.txt"
file_path4 = parent + "/results_GC/M31_H12/obj1/spec.txt"
file_path5 = parent + "/results_GC/M31_PANDAS_41/obj1/spec.txt"
file_path6 = parent + "/results_GC/Sextans_A_GC1/obj1/spec.txt"
files = [file_path1, file_path2, file_path3, file_path4, file_path5, file_path6]
names = ["DDO190", "F8D1", "M31 B336", "M31 H12", "M31 PA-41", "SextansA"]
redshifts = [0.000535, -0.000359, -0.0019824, -0.001262, -0.000351, 0.0009887]

n_files = len(files)
fig, axes = plt.subplots(n_files, 1, figsize=(7, 5*n_files), sharex=True)

for file, ax, name, z in zip(files, axes, names, redshifts):

    # Load the data
    data = np.loadtxt(file)
    x = data[:, 0]
    galaxy = data[:, 1]
    fit = data[:, 2]

    # Plot the data
    ax.plot(x, galaxy, color="black", drawstyle="steps-mid")
    ax.plot(x, fit, color="magenta", drawstyle="steps-mid")
    ax.margins(x=0, y=0.3)
    ax.set_title(name, fontsize=12, y=0.875, va='top')
    ax.tick_params(axis='both', which='major', labelsize=12)

    ax.axvline(x=4861.34 * (1 + z), color="grey", linestyle="dashed", linewidth=0.5)
    ax.axvline(x=4341.68 * (1 + z), color="grey", linestyle="dashed", linewidth=0.5)
    ax.axvline(x=4102.89 * (1 + z), color="grey", linestyle="dashed", linewidth=0.5)

    # dz = 4861.34 * (1 + z)
    # ax.set_xlim(dz*0.996, dz*1.004)
    
axes[0].text(x=(4861.34-50) * (1 + z), y=0, s="H$\u03B2$", fontsize=12, rotation='vertical', color="black")
axes[0].text(x=(4341.68-50) * (1 + z), y=0, s="H\u03B3", fontsize=12, rotation='vertical', color="black")
axes[0].text(x=(4102.89-50) * (1 + z), y=0, s="H\u03B4", fontsize=12, rotation='vertical', color="black")

xlim, xmax = axes[1].get_xlim()
ylim, ymax = axes[1].get_ylim()
for ax, file, z in zip(axes, files, redshifts):
    data = np.loadtxt(file)
    x = data[:, 0] / (1 + z)
    rect_start = patches.Rectangle((0, ylim), width=x[0], height=ymax*2, color='gray', alpha=0.35)
    ax.add_patch(rect_start)
    rect_end = patches.Rectangle((x[-1], ylim), width=xmax-x[-1], height=ymax*2, color='gray', alpha=0.35)
    ax.add_patch(rect_end)

plt.subplots_adjust(hspace=0)
fig.text(0.04, 0.5, 'Relative Flux', va='center', rotation='vertical', fontsize=12)
plt.xlabel("Wavelength $\AA$", fontsize=12)
plt.show()