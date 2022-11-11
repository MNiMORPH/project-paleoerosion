# Cosmogenic radionuclide concentrations

import pandas as pd
import numpy as np

# Input file
# Columns:
# CHECK EROSION RATE DEFINITON; NOW ENDING AT THIS TIME
# Age, Erosion rate ending at this time and starting from the last time
# Age, Erosion rate starting at this time and continuing until the next time
data = pd.read_csv('test_time_series.csv')

# Variables
P0 = 5. # atoms/g/yr: surface production rate
Lambda = 0.8 # meters: attenuation length

dt = data['Age [yr BP]'].diff(periods=-1)
#dt = -data['Age [yr BP]'].diff()

# Relative elevation of surface with respect to initial time
dz_surf = -dt * data['Erosion rate [mm/yr]']/1E3
dz_surf = np.roll(dz_surf, 1)
dz_surf[0] = 0
data['Relative surface elevation [m]'] = dz_surf.cumsum()

# New dt definition
#dt = -data['Age [yr BP]'].diff(periods=1)

# CRN concentrations
data['10Be concentration at surface [atoms/g]'] = pd.Series(dtype='float')

for i in data.index: #range(1, len(data.index)):
    # Depths, times, and erosion rates while sample is not yet eroded
    z_sample = data['Relative surface elevation [m]'][i] - \
                  data['Relative surface elevation [m]']
    _valid = z_sample < 0 # Only producing CRN when not yet eroded

    # Depths
    z_sample = z_sample[_valid]
    eros_rate = data['Erosion rate [mm/yr]'][_valid] / 1E3 # [m/yr]
    _dt = dt[_valid]

    # Piecewise linear integration
    dC = P0 * Lambda / eros_rate * ( 
                          np.exp( (z_sample + eros_rate * _dt) / Lambda ) - 
                          np.exp(z_sample / Lambda)
                          )

    # Summation
    data.loc[i, '10Be concentration at surface [atoms/g]'] = dC.sum()


### Plotting scratch space ###

from matplotlib import pyplot as plt
plt.plot(data['Age [yr BP]'][1:]/1000, data['10Be concentration at surface [atoms/g]'][1:], 'k-', linewidth=2)
plt.xlabel('Age [ka]', fontsize=14)
plt.ylabel('10Be concentration at surface [atoms/g]', fontsize=14)
plt.twinx()
plt.step(data['Age [yr BP]'][1:]/1000, data['Erosion rate [mm/yr]'][:-1], '0.5', linewidth=2)
plt.ylabel('Catchment-averaged erosion rate [mm/yr]', color='.5', fontsize=14)
plt.tight_layout()
plt.show()

