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

ages = np.array( _model_io['Age [yr BP]'] )

erosion_rate_min = 0.01 * np.ones(len(ages) - 1)
erosion_rate_max = 0.2 * np.ones(len(ages) - 1)

# Last time step constrained by data
erosion_rate_min[-1] = 0.3
erosion_rate_max[-1] = 0.5

cm = CosmicMonty(5, 0.7, ages, 'test_synth_data.csv')

cm.initialize_minmax_mode(erosion_rate_min, erosion_rate_max)
cm.initialize_ages(ages)
cm.initialize_output( csv_dir_2sigma = 
                              '/home/awickert/Desktop/CosmicMonty2sigmaTest' )

cm.mcloop(2)



