#! /usr/bin/python3

import sys
sys.path.append('./cosmicErosion')
import numpy as np
from cosmic_erosion import CosmicMCMC

# Fixed time grid [yr BP], oldest first.
# Starts at 100 ka to spin up the 10Be inventory before the oldest datum
# (~74 ka). The memory time is ~Λ/ε ≈ 0.7/0.0003 ≈ 2300 yr, so 26 ka of
# spin-up (100→74 ka) is well beyond saturation.
ages = np.array([100000, 80000, 70000, 60000, 50000, 40000,
                  30000, 25000, 20000, 15000, 10000,  5000, 0])

cm = CosmicMCMC(
    P0                = 17.5,
    attenuation_length= 0.7,
    ages              = ages,
    crn_data          = 'Mariotti2021.csv',
    erosion_rate_min  = 0.001,
    erosion_rate_max  = 2.0,
)

sampler = cm.run(n_steps=2000, progress=True)

# Diagnostics: mean acceptance fraction (healthy range ~0.2–0.5)
print('Mean acceptance fraction: {:.3f}'.format(
    np.mean(sampler.acceptance_fraction)))

# Save posterior samples and summary statistics before any plotting,
# so data is preserved even if the display is closed forcibly.
cm.save_chain('Mariotti_MCMC_chain.csv', discard=500, thin=10)
cm.save_summary_stats('Mariotti_MCMC_summary.csv', discard=500, thin=10)
print('Chain saved to Mariotti_MCMC_chain.csv')
print('Summary statistics saved to Mariotti_MCMC_summary.csv')

# Combined two-panel figure: pre-computes all model curves before
# opening the window, then shows once and closes cleanly.
cm.plot_summary(discard=500, thin=10, n_draws=300,
                savepath='Mariotti_MCMC_summary.png')
print('Figure saved to Mariotti_MCMC_summary.png')
