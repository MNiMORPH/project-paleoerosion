# Cosmogenic radionuclide concentrations

import pandas as pd
import numpy as np


class CosmicErosion(object):

    def __init__(self, data, P0, attenuation_length ):
        """
        :param [data]: Input data table with two columns named:
                        (1) Age [yr BP]
                        (2) Erosion rate [mm/yr]
                        The erosion rate is that starting at the associated
                        age and continuing until the subsequent time step.
                        This may be a Pandas DataFrame or a path to a CSV file.
        :type [data]: pandas.core.frame.DataFrame or csv
        
        :param [P0]: Surface 10Be production rate
        :type [P0]: float
        
        :param [attenuation_length]: 10Be production attenuation length [m]
        :type [attenuation_length]: float
        """

        if type(data) is str:
            self.data = pd.read_csv(data)
        elif type(data) is pandas.core.frame.DataFrame:
            self.data = data
        else:
            raise TypeError('Type must be str or Pandas DataFrame')
        
        self.P0 = P0
        self.attenuation_length = attenuation_length
        
    def initialize(self):
        """
        Set up the variables that are defined only once
        """

        # Relative surface elevation with time
        self.dt = self.data['Age [yr BP]'].diff(periods=-1)
        dz_surf = -self.dt * self.data['Erosion rate [mm/yr]']/1E3
        dz_surf = np.roll(dz_surf, 1)
        dz_surf[0] = 0
        self.data['Relative surface elevation [m]'] = dz_surf.cumsum()
        
        # Create column for 10Be concentration
        self.data['Modeled surface [10Be] [atoms/g]'] = pd.Series(dtype='float')
        
    def update(self, i):
        """
        Compute surface 10Be concentration as a time step
        :param [i]: Time-step index
        :type [i]: int
        """
        # Depths, times, and erosion rates while sample is not yet eroded
        z_sample = self.data['Relative surface elevation [m]'][i] - \
                      self.data['Relative surface elevation [m]']
        _valid = z_sample < 0 # Only producing CRN when not yet eroded

        # Depths
        z_sample = z_sample[_valid]
        eros_rate = self.data['Erosion rate [mm/yr]'][_valid] / 1E3 # [m/yr]
        _dt = self.dt[_valid]

        # Piecewise linear integration
        dC = self.P0 * self.attenuation_length / eros_rate * ( 
                              np.exp( (z_sample + eros_rate * _dt) / 
                                                self.attenuation_length ) - 
                              np.exp(z_sample / self.attenuation_length)
                              )

        # Summation
        self.data.loc[i, 'Modeled surface [10Be] [atoms/g]'] = dC.sum()

    def run(self):
        """
        Solve the 10Be concentration (via update()) for all time steps
        """
        for i in self.data.index: #range(1, len(data.index)):
            self.update(i)

    def finalize(self, show_plot=False, plot_savepath=None, csv_savepath=None):
        """
        Generate plots and/or output file(s)
        """
        if show_plot or plot_savepath:
            self.plot(show_plot, plot_savepath)
        if csv_savepath:
            self.data.to_csv(csv_savepath)

    def plot(self, show=True, savepath=None):
        """
        Plot model output
        """
        from matplotlib import pyplot as plt
        plt.plot(self.data['Age [yr BP]'][1:]/1000, self.data['Modeled surface [10Be] [atoms/g]'][1:], 'k-', linewidth=2)
        plt.xlabel('Age [ka]', fontsize=14)
        plt.ylabel('Modeled surface [10Be] [atoms/g]', fontsize=14)
        plt.twinx()
        plt.step(self.data['Age [yr BP]'][1:]/1000, self.data['Erosion rate [mm/yr]'][:-1], '0.5', linewidth=2)
        plt.ylabel('Catchment-averaged erosion rate [mm/yr]', color='.5', fontsize=14)
        plt.tight_layout()
        if savepath:
            plt.savefig(savepath)
        if show:
            plt.show()
        else:
            plt.close()
    
