#! /usr/bin/python3

from os import path
import pandas as pd
import numpy as np

#csv_dir = '/home/awickert/Desktop/CosmicMonty2sigmaTest'
csv_dir = '/home/awickert/Desktop/CosmicMontyTest04/2_sigma'

crn_data = pd.read_csv('test_synth_data.csv')

import glob
runfiles = sorted(glob.glob(path.join(csv_dir, '*.csv')))

from matplotlib import pyplot as plt
plt.ion()

fig = plt.figure(figsize=(6,12))

ax1 = plt.subplot(3,1,1)
ax2 = plt.subplot(3,1,2)
ax3 = plt.subplot(3,1,3)

for runfile in runfiles:
    model_io = pd.read_csv(runfile)
    ax1.plot( model_io['Age [yr BP]']/1000,
              model_io['Modeled surface [10Be] [atoms/g]'],
              'k-', linewidth=0.5, alpha=0.2 )
    # Integral of erosion rate: Should be smoother than local rate
    # Could also represent in mean, SD
    _tmp = np.cumsum( model_io['Erosion rate [mm/yr]'] * \
                      model_io['Age [yr BP]'].diff(-1) )
    _tmp = np.roll(_tmp,1)
    _tmp[0] = 0
    ax2.plot( model_io['Age [yr BP]']/1000, _tmp-_tmp[-1],
              'k-', linewidth=0.5, alpha=0.2 )

if crn_data is not None:
    ax1.errorbar( crn_data['Age [yr BP]']/1000,
                  crn_data['10Be Concentration [atoms/g]'],
                  yerr=crn_data['10Be SD [atoms/g]'],
                  marker='o', color='k', linestyle='', elinewidth=2 )



# Inefficient but whatever
rflist = []
for runfile in runfiles:
    rflist.append( pd.read_csv(runfile) )

erosResults = pd.DataFrame(columns=['Age [yr BP]', 'mean', 'median', 'sd'])
erosResults['Age [yr BP]'] = model_io['Age [yr BP]']
for i in model_io.index:
    _tmp = []
    for df in rflist:
        _tmp.append(df['Erosion rate [mm/yr]'][i])
    erosResults.loc[i, 'mean'] = np.mean(_tmp)
    erosResults.loc[i, 'median'] = np.median(_tmp)
    erosResults.loc[i, 'sd'] = np.std(_tmp)

# Assume constant through present -- no further information
# For the purpose of plotting
erosResults.loc[40, 'mean'] = erosResults.loc[39, 'mean']
erosResults.loc[40, 'median'] = erosResults.loc[39, 'median']
erosResults.loc[40, 'sd'] = erosResults.loc[39, 'sd']


# Use standard error
ax3.step( erosResults['Age [yr BP]']/1000., erosResults['median'], 'k-',
          where='post')
ax3.fill_between( np.array(erosResults['Age [yr BP]']/1000., dtype=float),
                  np.array(erosResults['mean'] - 
                  erosResults['sd']/len(erosResults['sd'])**.5, dtype=float),
                  np.array(erosResults['mean'] + 
                  erosResults['sd']/len(erosResults['sd'])**.5, dtype=float),
                  step='post',
                  color='k',
                  alpha=0.5 )
# SD
#ax3.fill_between( np.array(erosResults['Age [yr BP]']/1000., dtype=float),
#                  np.array(erosResults['mean']-erosResults['sd'], dtype=float),
#                  np.array(erosResults['mean']+erosResults['sd'], dtype=float),
#                  step='post',
#                  color='k',
#                  alpha=0.5 )

ax1.set_xlim(( model_io['Age [yr BP]'][len(model_io)-1]/1000, model_io['Age [yr BP]'][1]/1000 + 5 ))
ax2.set_xlim(( model_io['Age [yr BP]'][len(model_io)-1]/1000, model_io['Age [yr BP]'][1]/1000 + 5 ))
ax3.set_xlim(( model_io['Age [yr BP]'][len(model_io)-1]/1000, model_io['Age [yr BP]'][1]/1000 + 5 ))

ax1.set_xlim(-.5,25)
ax1.set_ylim(50000,90000)
ax2.set_xlim(-.5,25)
ax2.set_ylim(-1500,0)
ax3.set_xlim(-.5,25)

ax1.set_ylabel('Modeled surface\n[10Be] [atoms/g]', fontsize=14)
ax2.set_ylabel('Denudation [mm]', fontsize=14)
ax3.set_ylabel('Catchment-averaged\nerosion rate [mm/yr]', fontsize=14)
ax3.set_xlabel('Age [ka]', fontsize=14)

plt.tight_layout()

plt.show()

