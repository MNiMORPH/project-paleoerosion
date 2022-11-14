# Cosmogenic radionuclide concentrations

import pandas as pd
import numpy as np


class CosmicErosion(object):
    """
    Simulate the cosmogenic 10Be concentration of sediments gathered from
    a drainage basin in order to relate it to a possible history of
    catchment-averaged erosion rate.
    """

    def __init__(self, model_io, P0, attenuation_length, crn_data = None ):
        """
        :param [model_io]: Input table with two columns named:
                        (1) Age [yr BP]
                        (2) Erosion rate [mm/yr]
                        The erosion rate is that starting at the associated
                        age and continuing until the subsequent time step.
                        This may be a Pandas DataFrame or a path to a CSV file.
        :type [model_io]: pandas.core.frame.DataFrame or csv
        
        :param [P0]: Surface 10Be production rate
        :type [P0]: float
        
        :param [attenuation_length]: 10Be production attenuation length [m]
        :type [attenuation_length]: float

        :param [crn_data]: Input table with three columns named:
                        (1) Age [yr BP]
                        (2) 10Be Concentration [atoms/g]
                        (3) 10Be SD [atoms/g]
                        where "SD" is the standard deviation.
        :type [crn_data]: pandas.core.frame.DataFrame or str or None, optional
        """

        if type(model_io) is str:
            self.model_io = pd.read_csv(model_io)
        elif type(model_io) is pd.core.frame.DataFrame:
            self.model_io = model_io
        else:
            raise TypeError('Type must be str or Pandas DataFrame')
        
        self.P0 = P0
        self.attenuation_length = attenuation_length
        
        # Check for CRN data
        if type(crn_data) is str:
            self.crn_data = pd.read_csv(crn_data)
        elif type(crn_data) is pd.core.frame.DataFrame:
            self.crn_data = crn_data
        elif crn_data is None:
            pass
        else:
            raise TypeError('Type must be str or Pandas DataFrame or None')
    
    #######################
    # GETTERS AND SETTERS #
    #######################

    def set_P0(self, P0):
        self.P0 = P0
        
    def get_P0(self):
        return self.P0

    def set_attenuation_length(self):
        self.attenuation_length = attenuation_length
        
    def get_attenuation_length(self):
        return self.attenuation_length

    def set_model_io(self):
        self.model_io = model_io
    
    def get_model_io(self):
        return self.model_io
        
    ##########################
    # IRUF interface (CSDMS) #
    ##########################

    def initialize(self):
        """
        Set up the variables that are defined only once
        """

        # Relative surface elevation with time
        self.dt = self.model_io['Age [yr BP]'].diff(periods=-1)
        dz_surf = -self.dt * self.model_io['Erosion rate [mm/yr]']/1E3
        dz_surf = np.roll(dz_surf, 1)
        dz_surf[0] = 0
        self.model_io['Relative surface elevation [m]'] = dz_surf.cumsum()
        
        # Create column for 10Be concentration
        self.model_io['Modeled surface [10Be] [atoms/g]'] = pd.Series(dtype='float')
        
    def update(self, i):
        """
        Compute surface 10Be concentration as a time step
        :param [i]: Time-step index
        :type [i]: int
        """
        # Depths, times, and erosion rates while sample is not yet eroded
        z_sample = self.model_io['Relative surface elevation [m]'][i] - \
                      self.model_io['Relative surface elevation [m]']
        _valid = z_sample < 0 # Only producing CRN when not yet eroded

        # Depths
        z_sample = z_sample[_valid]
        eros_rate = self.model_io['Erosion rate [mm/yr]'][_valid] / 1E3 # [m/yr]
        _dt = self.dt[_valid]

        # Piecewise linear integration
        dC = self.P0 * self.attenuation_length / eros_rate * ( 
                              np.exp( (z_sample + eros_rate * _dt) / 
                                                self.attenuation_length ) - 
                              np.exp(z_sample / self.attenuation_length)
                              )

        # Summation
        self.model_io.loc[i, 'Modeled surface [10Be] [atoms/g]'] = dC.sum()

    def run(self):
        """
        Solve the 10Be concentration (via update()) for all time steps
        """
        for i in self.model_io.index:
            self.update(i)

    def finalize(self, show_plot=False, plot_savepath=None, csv_savepath=None):
        """
        Generate plots and/or output file(s)
        """
        if show_plot or plot_savepath:
            self.plot(show_plot, plot_savepath)
        if csv_savepath:
            self.model_io.to_csv(csv_savepath)

    ############
    # PLOTTING #
    ############
    
    def plot(self, show=True, savepath=None):
        """
        Plot model output
        """
        from matplotlib import pyplot as plt
        plt.plot(self.model_io['Age [yr BP]'][1:]/1000, self.model_io['Modeled surface [10Be] [atoms/g]'][1:], 'k-', linewidth=2)
        plt.xlabel('Age [ka]', fontsize=14)
        plt.ylabel('Modeled surface [10Be] [atoms/g]', fontsize=14)
        plt.twinx()
        plt.step(self.model_io['Age [yr BP]'][1:]/1000, self.model_io['Erosion rate [mm/yr]'][:-1], '0.5', linewidth=2)
        plt.ylabel('Catchment-averaged erosion rate [mm/yr]', color='.5', fontsize=14)
        plt.tight_layout()
        if savepath:
            plt.savefig(savepath)
        if show:
            plt.show()
        else:
            plt.close()
    
