#! /usr/bin/python3

import sys
sys.path.append('./cosmicErosion')
import numpy as np
from cosmic_erosion import CosmicMCMC

# 3000-year time grid [yr BP], oldest first.
# 28 ages → 27 free erosion rate parameters.
# Starts at 80 ka (explicit, since 80000 is not a multiple of 3000);
# remaining steps are uniform 3-kyr intervals down to 0.
# The steady-state spin-up correction (n_oldest=1) sets the initial 10Be
# inventory from the oldest datum so no long spin-up needed.
ages = np.hstack([[80000], np.arange(0, 78001, 3000)[::-1]])

cm = CosmicMCMC(
    P0                = 17.5,
    attenuation_length= 0.7,
    ages              = ages,
    crn_data          = 'Mariotti2021.csv',
    erosion_rate_min  = 0.001,
    erosion_rate_max  = 2.0,
    n_oldest          = 1,
    # Smoothness prior: sigma=0.5 → ~1.6x change per step at 1σ.
    # Suppresses unconstrained oscillations between data-sparse intervals.
    # IMPORTANT: high→low erosion transitions appear more readily in the
    # posterior than low→high. A drop in rate adds 10Be near the surface
    # quickly; a rise in rate must first erode away the entire pre-existing
    # subsurface inventory before the concentration record responds. Do not
    # over-interpret apparent rate increases.
    sigma_log_rate    = 0.5,
)

print('Parameters: {:d}   Walkers: {:d}'.format(
    cm.n_params, max(4 * cm.n_params, 16)))


sampler = cm.run(n_steps=2000, progress=True)

print('Mean acceptance fraction: {:.3f}'.format(
    np.mean(sampler.acceptance_fraction)))

# _3ka suffix avoids overwriting existing run outputs.
cm.save_chain('Mariotti_MCMC_3ka_80ka_chain.csv', discard=500, thin=10)
cm.save_summary_stats('Mariotti_MCMC_3ka_80ka_summary.csv', discard=500, thin=10)
print('Chain saved to Mariotti_MCMC_3ka_80ka_chain.csv')
print('Summary statistics saved to Mariotti_MCMC_3ka_80ka_summary.csv')

cm.plot_summary(discard=500, thin=10, n_draws=300,
                savepath='Mariotti_MCMC_3ka_80ka_summary.png')
print('Figure saved to Mariotti_MCMC_3ka_80ka_summary.png')
