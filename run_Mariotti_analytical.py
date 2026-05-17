#! /usr/bin/python3

import sys
sys.path.append('./cosmicErosion')
import numpy as np
from cosmic_erosion import CosmicAnalytical

ca = CosmicAnalytical(
    P0               = 17.5,
    attenuation_length = 0.7,
    crn_data         = 'Mariotti2021.csv',
    erosion_rate_min = 0.001,
    erosion_rate_max = 2.0,
    # Samples within 200 yr of each other are merged into an
    # inverse-variance-weighted mean before solving. Only the
    # {19773, 19763, 19603} yr BP triplet requires merging: the
    # 10-yr gap between the first two makes separate nodes physically
    # impossible (max accumulation in 10 yr is P0*dt = 175 atoms/g,
    # far less than the observed 17,680-atom jump). All other pairs
    # are >290 yr apart and are left as individual nodes.
    cluster_dt       = 200,
    # cluster_dt       = 0,
)

if ca.cluster_dt:
    print('Clustered data points: {:d}  (from {:d} raw samples)'.format(
        len(ca.ages_clustered), len(ca.crn_data)))
else:
    print('Data points: {:d}  (clustering off; raw samples: {:d})'.format(
        len(ca.ages_clustered), len(ca.crn_data)))

# Point-estimate solve: one erosion rate per clustered datum.
rates, ages = ca.solve()
print('\nPoint-estimate erosion rates:')
print('  {:>12s}  {:>12s}'.format('Age [ka BP]', 'Rate [mm/yr]'))
for age, rate in zip(ages / 1000, rates):
    print('  {:>12.1f}  {:>12.3f}'.format(age, rate))

# Monte Carlo uncertainty propagation: perturb concentrations by their
# 1-sigma errors and repeat the sequential solve for each replicate.
print('\nPropagating uncertainty (n=5000)...')
ca.propagate_uncertainty(n_mc=5000)

ca.save_summary('Mariotti_analytical_summary.csv')
print('Summary statistics saved to Mariotti_analytical_summary.csv')

ca.plot_summary(savepath='Mariotti_analytical_summary.png')
print('Figure saved to Mariotti_analytical_summary.png')

ca.plot_comparison(savepath='Mariotti_analytical_comparison.png')
print('Figure saved to Mariotti_analytical_comparison.png')
