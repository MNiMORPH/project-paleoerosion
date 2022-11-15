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
_model_io = pd.read_csv('test_time_series.csv')

#ages = np.array( _model_io['Age [yr BP]'] )
ages = np.arange(0,22001,1000)
ages = np.hstack((ages, [140000]))
ages = ages[::-1]

ages = np.arange(0,40000,1000)
ages = np.hstack((ages, [140000]))
ages = ages[::-1]

erosion_rate_min = 0.005 * np.ones(len(ages) - 1)
erosion_rate_max = 0.1 * np.ones(len(ages) - 1)

# Last time step constrained by data
erosion_rate_min[-1] = 0.005
erosion_rate_max[-1] = 0.1

# Spin up should be close to initial concentration
# Maybe make a function to help wtih this
#erosion_rate_min[0] = 0.065
#erosion_rate_max[0] = 0.07

cm = CosmicMonty(5, 0.7, ages, 'test_synth_data.csv')

cm.initialize_minmax_mode(erosion_rate_min, erosion_rate_max)
cm.initialize_ages(ages)
cm.initialize_output( csv_dir_1sigma = 
                      '/home/awickert/Desktop/CosmicMontyTest04/1_sigma',
                        csv_dir_2sigma = 
                      '/home/awickert/Desktop/CosmicMontyTest04/2_sigma',
                     )

cm.mcloop(50000, verbose=True)



