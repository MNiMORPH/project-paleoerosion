#! /usr/bin/python3

import sys
sys.path.append('./cosmicErosion')
import numpy as np
import pandas as pd
from cosmic_erosion import CosmicAnalytical

ca = CosmicAnalytical(
    P0               = 5,
    attenuation_length = 0.8,
    crn_data         = 'test_synth_data.csv',
    erosion_rate_min = 0.001,
    erosion_rate_max = 2.0,
    cluster_dt       = 0,   # samples well-separated; no clustering needed
)

print('Data points: {:d}'.format(len(ca.ages_clustered)))

rates, ages = ca.solve()
print('\nPoint-estimate erosion rates:')
print('  {:>12s}  {:>12s}'.format('Age [ka BP]', 'Rate [mm/yr]'))
for age, rate in zip(ages / 1000, rates):
    print('  {:>12.1f}  {:>12.4f}'.format(age, rate))

# Load true time series for comparison
truth = pd.read_csv('test_time_series.csv').dropna()
print('\nTrue erosion rates (test_time_series.csv):')
print('  {:>12s}  {:>12s}'.format('Age [ka BP]', 'Rate [mm/yr]'))
for _, row in truth.iterrows():
    print('  {:>12.1f}  {:>12.4f}'.format(row['Age [yr BP]'] / 1000,
                                           row['Erosion rate [mm/yr]']))

print('\nPropagating uncertainty (n=5000)...')
ca.propagate_uncertainty(n_mc=5000)

ca.save_summary('synth_analytical_summary.csv')
print('Summary statistics saved to synth_analytical_summary.csv')

ca.plot_summary(savepath='synth_analytical_summary.png')
print('Figure saved to synth_analytical_summary.png')

ca.plot_comparison(savepath='synth_analytical_comparison.png')
print('Figure saved to synth_analytical_comparison.png')
