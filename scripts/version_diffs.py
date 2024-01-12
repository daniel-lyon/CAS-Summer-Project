import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-1100, 1100, 100)

# Velocity
vel_901 = [94, 62, 49, 85, 90, 100, 166, -65, -596, -355, -102, 297]
vel_6715 = [-75, -67, -66, -65, -68, -66, 175, -50, -18, -17, -70, -27]
vel_main = [156, 121, 106, 143, 154, 158, 179, -21, -579, -363, -32, 326]

# Velocity dispersion
sig_901 = [70, 1, 30, 29, 1, 16, 8, 19, 19, 28, 96, 13]
sig_6715 = [76, 1, 24, 29, 1, 1, 0, 26, 9, 11, 33, 6]
sig_main = [80, 1, 32, 24, 1, 15, 5, 1, 7, 9, 42, 0]

# Log age
age_901 = [9.98, 10, 9.79, 9.86, 9.89, 9.8, 9.16, 8.67, 10.1, 9.56, 10.1, 9.96]
age_6715 = [9.7, 9.9, 9.7, 10, 10.1, 9.87, 7.84, 9.65, 10.1, 10.1, 9.86, 9.77]

# Metallicity
mh_901 = [-1.4, -1.21, -1.71, -0.71, -0.624, -1.69, -1.6, -0.0772, -1.71, -1.64, -1.37, -1.71]
mh_6715 = [-1.44, -1.15, -1.92, -0.731, -0.655, -1.6, -2.27, -1.53, -2.19, -1.91, -1.43, -2.24]

# Create a figure with two subplots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 8))

# Velocity
ax1.plot(x, x, color='blue', label='y=x')
ax1.scatter(vel_6715, vel_901, color='red', label='v6.7.15 vs v9.0.1')
ax1.scatter(vel_main, vel_901, color='green', label='GC Fit vs SP Fit', marker='x')
ax1.set_xlabel('v6.7.15 Velocity (km/s)')
ax1.set_ylabel('v9.0.1 Velocity (km/s)')
ax1.set_xscale('symlog')
ax1.set_yscale('symlog')
ax1.set_xlim(-1e3, 1e3)
ax1.set_ylim(-1e3, 1e3)
ax1.grid()
ax1.legend()

# Sigma
ax2.plot(x, x, color='blue', label='y=x')
ax2.scatter(sig_6715, sig_901, color='red', label='v6.7.15 vs v9.0.1')
ax2.scatter(sig_main, sig_901, color='green', label='GC Fit vs SP Fit', marker='x')
ax2.set_xlabel('v6.7.15 Sigma (km/s)')
ax2.set_ylabel('v9.0.1 Sigma (km/s)')
ax2.set_xlim(-5, 100)
ax2.set_ylim(-5, 100)
ax2.grid()
ax2.legend()

# Log age
ax3.plot(x, x, color='blue', label='y=x')
ax3.scatter(age_6715, age_901, color='red', label='v6.7.15 vs v9.0.1')
ax3.set_xlabel('v6.7.15 Log Age')
ax3.set_ylabel('v9.0.1 Log Age')
ax3.set_xlim(7, 11)
ax3.set_ylim(7, 11)
ax3.grid()
ax3.legend()

# Metallicity
ax4.plot(x, x, color='blue', label='y=x')
ax4.scatter(mh_6715, mh_901, color='red', label='v6.7.15 vs v9.0.1')
ax4.set_xlabel('v6.7.15 Metallicity')
ax4.set_ylabel('v9.0.1 Metallicity')
ax4.set_xlim(-3, 1)
ax4.set_ylim(-3, 1)
ax4.grid()
ax4.legend()

# Adjust the spacing between subplots
plt.tight_layout()

# Show the plot
plt.show()
