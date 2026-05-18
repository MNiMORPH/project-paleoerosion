#! /usr/bin/python3

import sys
sys.path.append('./cosmicErosion')
import numpy as np
from cosmic_erosion import CosmicAnalytical

# Parameters from Penprase et al. (2025), Table S10.
# P0 = 6.41 atoms/g/yr (locally-scaled total production rate).
# Penprase et al. report SLHL spallogenic P0 = 5.38 at/g/yr and muogenic
# P0 = 0.196 at/g/yr (total SLHL = 5.576). Back-calculating from their
# published steady-state rates and Table S10 concentrations gives a
# consistent local scaling factor of 1.150, yielding P0_total_local =
# 6.41 at/g/yr. Their treatment folds both components into a single
# effective production rate with the spallogenic attenuation length
# (Λ = 0.714 m), so P0 = 6.41 with Λ = 0.714 exactly reproduces their
# published erosion rates.
# Λ = 0.714 m (locally-scaled effective attenuation length, from paper).

ca = CosmicAnalytical(
    P0               = 6.185,   # locally-scaled spallogenic only
    attenuation_length = 0.714,
    crn_data         = 'Penprase2025.csv',
    erosion_rate_min = 0.001,
    erosion_rate_max = 2.0,
    # Lower-Deep and Upper-Deep are 120 yr apart — merged by cluster_dt=200.
    cluster_dt       = 200,
    # Muogenic production (Penprase et al. 2025): SLHL 0.196 at/g/yr,
    # scaled by 1.150 to local conditions. Lambda_mu=5.6 m corresponds to
    # stopped muons (Heisinger et al. 2002, ~1500 g/cm² / 2.7 g/cm³).
    P0_mu            = 0.225,
    Lambda_mu        = 5.6,
)

if ca.cluster_dt:
    print('Clustered data points: {:d}  (from {:d} raw samples)'.format(
        len(ca.ages_clustered), len(ca.crn_data)))
else:
    print('Data points: {:d}  (clustering off; raw samples: {:d})'.format(
        len(ca.ages_clustered), len(ca.crn_data)))

rates, ages = ca.solve()
print('\nPoint-estimate erosion rates:')
print('  {:>12s}  {:>12s}'.format('Age [ka BP]', 'Rate [mm/yr]'))
for age, rate in zip(ages / 1000, rates):
    print('  {:>12.2f}  {:>12.4f}'.format(age, rate))

print('\nPropagating uncertainty (n=5000)...')
ca.propagate_uncertainty(n_mc=5000)

ca.save_summary('Penprase_analytical_summary.csv')
print('Summary statistics saved to Penprase_analytical_summary.csv')

ca.plot_summary(savepath='Penprase_analytical_summary.png')
print('Figure saved to Penprase_analytical_summary.png')

ca.plot_comparison(savepath='Penprase_analytical_comparison.png')
print('Figure saved to Penprase_analytical_comparison.png')
