#! /usr/bin/python3

import sys
sys.path.append('./cosmicErosion')
from cosmic_erosion import CosmicErosion

ce = CosmicErosion('test_time_series.csv', 5, 0.7,
                    crn_data='test_synth_data.csv')
ce.initialize()
ce.run()
ce.evaluate()
ce.finalize(show_plot=True)
#ce.finalize(show_plot=False, plot_savepath='test2.pdf', csv_savepath='test2.csv')

"""
# Scratch for future
ce.evaluate() # compare with data
if (ce.crn['within2SDbla'] == True).all()
    ce.finalize(show_plot=False, plot_savepath='test2.pdf', csv_savepath='test2.csv')
else:
    ce.finalize(show_plot=False)
"""
