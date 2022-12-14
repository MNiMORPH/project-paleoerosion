# Cosmogenic radionuclide concentrations

import pandas as pd
import numpy as np
from os import path
import warnings
from scipy.stats import gamma, beta

#################
# FORWARD MODEL #
#################

class CosmicErosion(object):
    """
    Simulate the cosmogenic 10Be concentration of sediments gathered from
    a drainage basin in order to relate it to a possible history of
    catchment-averaged erosion rate.
    """
    import pandas as pd
    import numpy as np
    from os import path

    def __init__(self, model_io, P0, attenuation_length, crn_data = None ):
        """
        :param [model_io]: Input table with two columns named:
                        (1) Age [yr BP]
                        (2) Erosion rate [mm/yr]
                        The erosion rate is that starting at the associated
                        age and continuing until the subsequent time step.
                        This may be a Pandas DataFrame or a path to a CSV file.
                        It can be set to None, in which case it will be created
                        as an empty DataFrame that can be populated later.
        :type [model_io]: pandas.core.frame.DataFrame or str or None
        
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
        elif type(model_io) is type(None):
            self.model_io = pd.DataFrame( columns=['Age [yr BP]', 
                                                   'Erosion rate [mm/yr]'] )
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
            self.crn_data = crn_data
        else:
            raise TypeError('Type must be str or Pandas DataFrame or None')
    
    #######################
    # GETTERS AND SETTERS #
    #######################

    def set_P0(self, P0):
        self.P0 = P0
        
    def get_P0(self):
        return self.P0

    def set_attenuation_length(self, attenuation_length):
        self.attenuation_length = attenuation_length
        
    def get_attenuation_length(self):
        return self.attenuation_length

    def set_model_io(self, model_io):
        self.model_io = model_io
    
    def get_model_io(self):
        return self.model_io
        
    def set_crn_data(self, crn_data):
        self.crn_data = crn_data
    
    def get_crn_data(self):
        return self.crn_data

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
                              np.exp( np.array(z_sample + eros_rate * _dt,
                                                dtype=float) / 
                                                self.attenuation_length ) - 
                              np.exp( np.array(z_sample, dtype=float) / 
                                                self.attenuation_length )
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

    ###############################
    # EVALUATE MODEL AGAINST DATA #
    ###############################
    
    def evaluate(self):
        """
        Evaluate model against data
        """
        
        if self.crn_data is None:
            raise TypeError('None-type crn_data: Should be Pandas DataFrame')
        
        #diff = crn.merge(data[['Age [yr BP]','10Be concentration at surface [atoms/g]']], on='Age [yr BP]')
        
        # In case join is faster
        self.crn_data = self.crn_data.join( 
                      self.model_io.set_index('Age [yr BP]')
                      ['Modeled surface [10Be] [atoms/g]'], on='Age [yr BP]'
                      )
        self.crn_data['10Be error [atoms/g]'] = np.abs( 
                      self.crn_data['10Be Concentration [atoms/g]'] - 
                      self.crn_data['Modeled surface [10Be] [atoms/g]']
                      )
        self.crn_data['Within 2SD'] = self.crn_data['10Be error [atoms/g]'] < \
                                            2*self.crn_data['10Be SD [atoms/g]']
        self.crn_data['Within 1SD'] = self.crn_data['10Be error [atoms/g]'] < \
                                              self.crn_data['10Be SD [atoms/g]']

    ############
    # PLOTTING #
    ############
    
    def plot(self, show=True, savepath=None):
        """
        Plot model output
        """
        from matplotlib import pyplot as plt
        plt.plot(self.model_io['Age [yr BP]'][1:]/1000, self.model_io['Modeled surface [10Be] [atoms/g]'][1:], 'k-', linewidth=2)
        if self.crn_data is not None:
            plt.errorbar( self.crn_data['Age [yr BP]']/1000,
                          self.crn_data['10Be Concentration [atoms/g]'],
                          yerr=self.crn_data['10Be SD [atoms/g]'],
                          marker='o', color='k', linestyle='', elinewidth=2 )
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
    

################################
# MONTE CARLO INVERSE MODELING #
################################

class CosmicMonty(object):

    def __init__(self, P0, attenuation_length, ages, crn_data):
        # !!!!!!!! WARNING: AGES NOT CURRENTLY USED HERE <-- FIX LATER
        self.SD_mode = False
        self.MinMax_mode = False
        self.ce = CosmicErosion(None, P0, attenuation_length, crn_data=crn_data)
        # Placeholder for including the gamma-function scaling parameter
        # scale=?
        # shape=1/scale [* mean of last data point]
        
    def initialize_minmax_mode(self, erosion_rate_min, erosion_rate_max):
        self.SD_mode = False
        self.MinMax_mode = True
        self.erosion_rate_min = erosion_rate_min
        self.erosion_rate_max = erosion_rate_max

    def initialize_SD_mode(self, erosion_rate_mean, erosion_rate_SD):
        self.SD_mode = False
        self.MinMax_mode = True
        self.erosion_rate_mean = erosion_rate_mean
        self.erosion_rate_SD = erosion_rate_SD

    def initialize_ages(self, ages):
        # Reaching directly into table instead of using any kind of setter
        self.ce.model_io['Age [yr BP]'] = ages

    def initialize_output(self, csv_dir_1sigma=None, plot_dir_1sigma=None,
                                csv_dir_2sigma=None, plot_dir_2sigma=None ):
        self.csv_dir_1sigma = csv_dir_1sigma
        self.plot_dir_1sigma = plot_dir_1sigma
        self.csv_dir_2sigma = csv_dir_2sigma
        self.plot_dir_2sigma = plot_dir_2sigma

    def mcloop(self, n, verbose=False):
        for i in range(n):
            npad = len(str(n))
            if verbose:
                print( ('%'+str(npad)+'s') %(i+1) + ' / ' + str(n) + ' -- ',
                        end='' )
            if self.SD_mode:
                # Should work with scalars or Numpy arrays
                self.ce.model_io.loc[ :len(self.ce.model_io) - 2, 
                                        'Erosion rate [mm/yr]' ] = \
                    self.erosion_rate_SD \
                    * np.random.standard_normal( len(self.ce.model_io
                                                      ['Age [yr BP]']) - 1 ) \
                    + self.erosion_rate_mean
            elif self.MinMax_mode:
                # Should work with scalars or Numpy arrays
                self.ce.model_io.loc[ :len(self.ce.model_io) - 2, 
                                        'Erosion rate [mm/yr]' ] = \
                    (self.erosion_rate_max - self.erosion_rate_min) \
                    * np.random.random_sample( len(self.ce.model_io
                                                      ['Age [yr BP]']) - 1 ) \
                    + self.erosion_rate_min
            else:
                raise ValueError("Neither SD_mode nor MinMax_mode are set.")
            
            self.ce.initialize()
            self.ce.run()
            self.ce.evaluate()
            if verbose:
                print( np.sum((self.ce.crn_data['Within 2SD'] == True)),
                       'of', len(self.ce.crn_data),
                       'data points within error of model outputs.' )
            if (self.ce.crn_data['Within 1SD'] == True).all():
                if self.csv_dir_1sigma is not None:
                    self.ce.model_io.to_csv( path.join(self.csv_dir_1sigma,
                                'model_run_'+('%0'+str(npad)+'d') %i+'.csv') )
                if self.plot_dir_1sigma is not None:
                    self.ce.plot( show=False,
                                  savepath=path.join(self.csv_dir_1sigma,
                                'model_run_'+('%0'+str(npad)+'d') %i)+'.png' )
            if (self.ce.crn_data['Within 2SD'] == True).all():
                if self.csv_dir_2sigma is not None:
                    self.ce.model_io.to_csv( path.join(self.csv_dir_2sigma,
                                'model_run_'+('%0'+str(npad)+'d') %i+'.csv') )
                if self.plot_dir_2sigma is not None:
                    self.ce.plot( show=False,
                                  savepath=path.join(self.csv_dir_2sigma,
                                'model_run_'+('%0'+str(npad)+'d') %i+'.png') )
            # Clean columns
            self.ce.crn_data = self.ce.crn_data.drop(
                  columns = [ 'Modeled surface [10Be] [atoms/g]',
                              '10Be error [atoms/g]',
                              'Within 2SD', 'Within 1SD'] )

class MontyPlot(object):

    def __init__(self, csv_dir):
        import glob
        runfiles = glob.glob(path.join(csv_dir, '*.csv'))
    

#############################
# BAYESIAN INVERSE MODELING #
#############################

class CosmicBay(object):

    from scipy.stats import gamma, beta

    def __init__(self, P0, attenuation_length, ages, crn_data, memory=1):
        self.ce = CosmicErosion(None, P0, attenuation_length, crn_data=crn_data)
        self.ages=ages
        self.bayesian_memory=memory # How many time steps to average over
                                    # in Bayesian series
        if memory != 1:
            warnings.warn("WARNING! BAYESIAN MEMORY NOT YET INCORPORATED!")
        
        # Hard-code memory for now
        self.bayesian_memory = beta.pdf(np.linspace(0, 1, 11), 4, 1)
        self.bayesian_memory /= np.sum(self.bayesian_memory)
    
    def initialize_ages(self, ages=None):
        # Reaching directly into table instead of using any kind of setter
        if ages is None:
            # And if self.ages is none, well, whoops. -- not possible now
            self.ce.model_io['Age [yr BP]'] = self.ages
        else:
            self.ce.model_io['Age [yr BP]'] = ages
    
    def initialize_erosion_rates(self, spinup_min=None, spinup_max=None,
                                       special_min=None, special_max=None,
                                       special_i=None):
        """
        Bayesian everywhere except where it is random (spinup, special times)
        special_min and special_max should be either of length len(special_i)
        or length 1 (untested, but I think)
        """
        
        if spinup_min is None:
            spinup_min = self.spinup_min
        else:
            self.spinup_min = spinup_min
        if spinup_max is None:
            spinup_max = self.spinup_max
        else:
            self.spinup_max = spinup_max
        if special_min is None:
            special_min = self.special_min
        else:
            self.special_min = special_min
        if special_max is None:
            special_max = self.special_max
        else:
            self.special_max = special_max
        if special_i is None:
            special_i = self.special_i
        else:
            self.special_i = special_i

        erosion_rates = np.full(len(self.ages), np.nan)

        # For copy/paste tests
        #erosion_rates *= np.nan

        # Initial erosion rate
        erosion_rates[0] = (spinup_max - spinup_min) \
                                * np.random.random_sample() \
                                + spinup_min

        # Set erosion rates at special times
        # To allow array operations if needed
        if not np.isscalar(special_max):
            special_max = np.array(special_max)
        if not np.isscalar(special_min):
            special_min = np.array(special_min)
        if not np.isscalar(special_i):
            n_special = len(special_i)
            special_i = np.array(special_i)
        else:
            n_special = 1
        special_erosion_rates = (special_max - special_min) \
                                * np.random.random_sample(n_special) \
                                + special_min
        erosion_rates[special_i] = special_erosion_rates


        # Bayesian series
        bayesian_memory = self.bayesian_memory
        # -1 because the last applies to the future
        for i in range(len(erosion_rates)-1):
            if np.isnan(erosion_rates[i]):
                # Gamma function, always positive,
                # with mean at last value
                if i >= 11:
                    memory_mean = np.sum( bayesian_memory
                                            * erosion_rates[i-11:i] )
                else:
                    memory_mean = np.sum( bayesian_memory[-i:]
                                            * erosion_rates[:i] ) \
                                    / np.sum(bayesian_memory[-i:])
                                # CHECK FOR ANY BUG HERE: CODING QUICKLY
                    
                erosion_rates[i] = memory_mean*gamma.rvs(4)/4.

        """
        # Bayesian series
        # -1 because the last applies to the future
        for i in range(len(erosion_rates) - 1):
            if np.isnan(erosion_rates[i]):
                # Gamma function, always positive,
                # with mean at last value
                try:
                    erosion_rates[i] = erosion_rates[i-1]*gamma.rvs(4)/4.
                except:
                    erosion_rates[i] = 0 #gamma.rvs( 1E-8 )
                    warnings.warn("0 wall")
                    # May have crashed into a 0 wall
                    # Let it stay there for now
        """


        self.erosion_rates = erosion_rates

        
    def initialize_output(self, csv_dir_1sigma=None, plot_dir_1sigma=None,
                                csv_dir_2sigma=None, plot_dir_2sigma=None ):
        self.csv_dir_1sigma = csv_dir_1sigma
        self.plot_dir_1sigma = plot_dir_1sigma
        self.csv_dir_2sigma = csv_dir_2sigma
        self.plot_dir_2sigma = plot_dir_2sigma

    def loop(self, n, verbose=False):
        for i in range(n):
            npad = len(str(n))
            if verbose:
                print( ('%'+str(npad)+'s') %(i+1) + ' / ' + str(n) + ' -- ',
                        end='' )
                        
            # Will use these values instead of initial initialization
            self.initialize_erosion_rates()
            #print(self.erosion_rates)
                        
            self.ce.model_io[ 'Erosion rate [mm/yr]' ] = self.erosion_rates
            
            self.ce.initialize()
            self.ce.run()
            self.ce.evaluate()
            if verbose:
                print( np.sum((self.ce.crn_data['Within 2SD'] == True)),
                       'of', len(self.ce.crn_data),
                       'data points within error of model outputs.' )
            if (self.ce.crn_data['Within 1SD'] == True).all():
                if self.csv_dir_1sigma is not None:
                    self.ce.model_io.to_csv( path.join(self.csv_dir_1sigma,
                                'model_run_'+('%0'+str(npad)+'d') %i+'.csv') )
                if self.plot_dir_1sigma is not None:
                    self.ce.plot( show=False,
                                  savepath=path.join(self.csv_dir_1sigma,
                                'model_run_'+('%0'+str(npad)+'d') %i)+'.png' )
            if (self.ce.crn_data['Within 2SD'] == True).all():
                if self.csv_dir_2sigma is not None:
                    self.ce.model_io.to_csv( path.join(self.csv_dir_2sigma,
                                'model_run_'+('%0'+str(npad)+'d') %i+'.csv') )
                if self.plot_dir_2sigma is not None:
                    self.ce.plot( show=False,
                                  savepath=path.join(self.csv_dir_2sigma,
                                'model_run_'+('%0'+str(npad)+'d') %i+'.png') )
            # Clean columns
            self.ce.crn_data = self.ce.crn_data.drop(
                  columns = [ 'Modeled surface [10Be] [atoms/g]',
                              '10Be error [atoms/g]',
                              'Within 2SD', 'Within 1SD'] )
            
