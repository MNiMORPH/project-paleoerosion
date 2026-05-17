#! /usr/bin/python3

import sys
sys.path.append('./cosmicErosion')
import pandas as pd
import numpy as np
from cosmic_erosion import CosmicErosion, CosmicMonty

################################
# CONSTRUCT ARRAY IN SAME CODE #
################################

# Bring in to set a range of ages
_model_io = pd.read_csv('Mariotti2021.csv')

ages = np.array( _model_io['Age [yr BP]'] )
#ages = np.arange(0,80000,5000)
ages = np.hstack((ages, [80000]))
ages = ages[::-1]

erosion_rate_min = 0.15 * np.ones(len(ages) - 1)
erosion_rate_max = .5 * np.ones(len(ages) - 1)

# Spin up should be close to initial concentration
# Maybe make a function to help wtih this
erosion_rate_min[0] = 0.15
erosion_rate_max[0] = 0.35

cm = CosmicMonty(17.5, 0.7, ages, 'Mariotti2021.csv')

cm.initialize_minmax_mode(erosion_rate_min, erosion_rate_max)
cm.initialize_ages(ages)
cm.initialize_output( csv_dir_1sigma = 
                      '/home/awickert/Desktop/CosmicMontyTest04/1_sigma',
                        csv_dir_2sigma = 
                      '/home/awickert/Desktop/CosmicMontyTest04/2_sigma',
                     )

cm.mcloop(1000, verbose=True)



