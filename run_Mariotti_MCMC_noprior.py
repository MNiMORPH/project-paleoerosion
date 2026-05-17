#! /usr/bin/python3

import sys
sys.path.append('./cosmicErosion')
import numpy as np
from cosmic_erosion import CosmicMCMC

# Same 3-kyr / 80-ka grid as run_Mariotti_MCMC_1ka.py, but with no
# smoothness prior (sigma_log_rate=None). The MCMC is free to change
# rate between any adjacent intervals without penalty, making it
# directly comparable to the analytical sequential solver.
ages = np.hstack([[80000], np.arange(0, 78001, 3000)[::-1]])

cm = CosmicMCMC(
    P0                = 17.5,
    attenuation_length= 0.7,
    ages              = ages,
    crn_data          = 'Mariotti2021.csv',
    erosion_rate_min  = 0.001,
    erosion_rate_max  = 2.0,
    n_oldest          = 1,
    sigma_log_rate    = None,
)

print('Parameters: {:d}   Walkers: {:d}'.format(
    cm.n_params, max(4 * cm.n_params, 16)))

sampler = cm.run(n_steps=2000, progress=True)

print('Mean acceptance fraction: {:.3f}'.format(
    np.mean(sampler.acceptance_fraction)))

cm.save_chain('Mariotti_MCMC_noprior_chain.csv', discard=500, thin=10)
cm.save_summary_stats('Mariotti_MCMC_noprior_summary.csv', discard=500, thin=10)
print('Chain saved to Mariotti_MCMC_noprior_chain.csv')
print('Summary statistics saved to Mariotti_MCMC_noprior_summary.csv')

cm.plot_summary(discard=500, thin=10, n_draws=300,
                savepath='Mariotti_MCMC_noprior_summary.png')
print('Figure saved to Mariotti_MCMC_noprior_summary.png')
