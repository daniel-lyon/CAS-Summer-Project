import pandas as pd
import matplotlib.pyplot as plt

# Setup
file_path = "/home/daniel/Documents/Swinburne/ultra-diffuse-galaxies/results_GC/results.ods"
sheet_name = 2  # Sheet index starts from 0
columns = ['Blue Vel', 'Red Vel', 'All Vel']

# Read excel file
df = pd.read_excel(file_path, sheet_name=sheet_name, usecols=columns)
df = df.dropna(how='all')
df['Red Diff'] = df['Red Vel'] - df['All Vel']
df['Blue Diff'] = df['Blue Vel'] - df['All Vel']

# Axis values
blue_vel = df['Blue Vel']
red_vel = df['Red Vel']
all_vel = df['All Vel']
names = ['NGC_247 GCs obj1', 'NGC_247 GCs obj2', 'NGC_247 GCs obj3', 
         'NGC_247 GCs2 obj1', 'NGC_247 GCs2 obj2', 'NGC_247 GCs2 obj3',
         'DDO190', 
         'F8D1', 
         'M31_B336', 'M31_H12', 'M31_PANDAS_41',
         'Sextans_A_GC1']

# Plotting line scatter
plt.scatter(names, blue_vel, label='Blue Vel', color='blue', marker='o', alpha=0.5)
plt.scatter(names, red_vel, label='Red Vel', color='red', marker='s', alpha=0.5)
plt.scatter(names, all_vel, label='All Vel', color='green', marker='^', alpha=0.5)
plt.xticks(rotation=90)
plt.xlabel('Source')
plt.ylabel('Velocity (km/s)')
plt.title('Velocity Comparison')
plt.legend(loc='lower left')
plt.grid(axis='y')
plt.tight_layout()
plt.show()

# Plot differences
plt.scatter(names, df['Blue Diff'], label='Blue Vel', color='blue', marker='o', alpha=0.5)
plt.scatter(names, df['Red Diff'], label='Red Vel', color='red', marker='s', alpha=0.5)
plt.xticks(rotation=90)
plt.xlabel('Source')
plt.ylabel('$\Delta$Velocity (km/s)')
plt.title('Velocity Comparison')
plt.legend(loc='lower right')
plt.grid(axis='y')
plt.tight_layout()
plt.show()