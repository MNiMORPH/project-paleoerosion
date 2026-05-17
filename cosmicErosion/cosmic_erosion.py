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


##########################
# MCMC INVERSE MODELING  #
##########################

class CosmicMCMC(object):
    """
    MCMC inverse modeling for cosmogenic 10Be erosion rate reconstruction.

    Uses a fixed user-specified time grid and the emcee ensemble sampler.
    The grid should extend before the oldest datum to allow spin-up of the
    10Be inventory toward steady state.

    Requires: emcee  (pip install emcee)
    """

    def __init__(self, P0, attenuation_length, ages, crn_data,
                 erosion_rate_min=0.01, erosion_rate_max=5.0, n_oldest=1,
                 sigma_log_rate=None):
        """
        :param ages: Fixed time grid, oldest first [yr BP].
        :param crn_data: Path to CSV or DataFrame with columns:
                         'Age [yr BP]', '10Be Concentration [atoms/g]',
                         '10Be SD [atoms/g]'
        :param erosion_rate_min: Lower bound [mm/yr]. Scalar or array of
                                 length len(ages)-1.
        :param erosion_rate_max: Upper bound [mm/yr]. Scalar or array of
                                 length len(ages)-1.
        :param n_oldest: Number of oldest data points to average for the
                         steady-state spin-up concentration. The mean
                         concentration is used as C_ss at the model start,
                         with exponential attenuation for grains that were
                         at depth at that time.
        :param sigma_log_rate: Standard deviation of the Gaussian random-walk
                               smoothness prior on log(erosion rate). Controls
                               how much the rate is allowed to change between
                               adjacent intervals: sigma=0.5 → ~1.6x change
                               per step at 1σ; sigma=1.0 → ~2.7x. Set to None
                               (default) to use a flat prior with no smoothing.
                               NOTE: this prior suppresses unconstrained jumps
                               between data-sparse intervals, but it introduces
                               an asymmetry — high→low erosion transitions are
                               detectable more readily than low→high, because a
                               drop in erosion quickly adds 10Be near the
                               surface, whereas a rise in erosion must first
                               strip away the entire pre-existing subsurface
                               inventory before the signal of the new rate
                               emerges. Keep this in mind when interpreting
                               inferred rate increases.
        """
        self.ages = np.array(ages)
        self.n_params = len(self.ages) - 1
        self.ce = CosmicErosion(None, P0, attenuation_length, crn_data=crn_data)
        self.ce.model_io['Age [yr BP]'] = self.ages
        self.erosion_rate_min = np.broadcast_to(
            np.atleast_1d(erosion_rate_min), self.n_params).copy()
        self.erosion_rate_max = np.broadcast_to(
            np.atleast_1d(erosion_rate_max), self.n_params).copy()
        # Steady-state surface [10Be] estimated from the n_oldest oldest data
        # points. Used to initialise grains with pre-existing inventory rather
        # than spinning up from zero.
        crn_sorted = self.ce.crn_data.sort_values('Age [yr BP]',
                                                   ascending=False)
        self._C_spinup = float(
            crn_sorted['10Be Concentration [atoms/g]'].iloc[:n_oldest].mean())
        self.sigma_log_rate = sigma_log_rate

    def _run_forward(self, erosion_rates):
        """
        Set erosion rates, run the forward model, and add the steady-state
        spin-up correction so that grains start with pre-existing 10Be rather
        than zero inventory.

        C_pre[i] = C_spinup * exp(z[i] / Lambda)

        where z[i] is the (negative) depth of sample i at the model start.
        Grains collected later were deeper at t_start, so they receive a
        smaller correction.
        """
        self.ce.model_io.loc[:self.n_params - 1,
                             'Erosion rate [mm/yr]'] = erosion_rates
        self.ce.initialize()
        self.ce.run()
        z = self.ce.model_io['Relative surface elevation [m]'].values.astype(float)
        self.ce.model_io['Modeled surface [10Be] [atoms/g]'] += (
            self._C_spinup * np.exp(z / self.ce.attenuation_length))

    def log_prior(self, erosion_rates):
        if (np.any(erosion_rates < self.erosion_rate_min) or
                np.any(erosion_rates > self.erosion_rate_max)):
            return -np.inf
        if self.sigma_log_rate is not None:
            diffs = np.diff(np.log(erosion_rates))
            return -0.5 * np.sum(diffs ** 2) / self.sigma_log_rate ** 2
        return 0.0

    def log_likelihood(self, erosion_rates):
        self._run_forward(erosion_rates)
        model_ages = self.ce.model_io['Age [yr BP]'].values
        model_conc = self.ce.model_io[
            'Modeled surface [10Be] [atoms/g]'].values
        if not np.all(np.isfinite(model_conc)):
            return -np.inf
        # Interpolate model output to data ages (np.interp needs increasing x)
        sort_idx = np.argsort(model_ages)
        conc_at_data = np.interp(
            self.ce.crn_data['Age [yr BP]'].values,
            model_ages[sort_idx],
            model_conc[sort_idx])
        obs = self.ce.crn_data['10Be Concentration [atoms/g]'].values
        sigma = self.ce.crn_data['10Be SD [atoms/g]'].values
        return -0.5 * np.sum(((obs - conc_at_data) / sigma) ** 2)

    def log_posterior(self, erosion_rates):
        lp = self.log_prior(erosion_rates)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_likelihood(erosion_rates)

    def run(self, n_steps=1000, n_walkers=None, starting_guess=None,
            progress=True):
        """
        Run the MCMC sampler.

        :param n_steps: Steps per walker.
        :param n_walkers: Number of walkers (default: 4 * n_params).
        :param starting_guess: Initial erosion rates [mm/yr], length n_params.
                               Default: midpoint of prior bounds.
        :param progress: Show tqdm progress bar.
        :return: emcee EnsembleSampler (also stored as self.sampler).
        """
        import emcee
        if n_walkers is None:
            n_walkers = max(4 * self.n_params, 16)
        if starting_guess is None:
            starting_guess = 0.5 * (self.erosion_rate_min + self.erosion_rate_max)
        starting_guess = np.asarray(starting_guess)
        rng = np.random.default_rng()
        p0 = starting_guess * (
            1.0 + 0.05 * rng.standard_normal((n_walkers, self.n_params)))
        p0 = np.clip(p0, self.erosion_rate_min, self.erosion_rate_max)
        self.sampler = emcee.EnsembleSampler(
            n_walkers, self.n_params, self.log_posterior)
        self.sampler.run_mcmc(p0, n_steps, progress=progress)
        return self.sampler

    def get_flat_chain(self, discard=100, thin=10):
        """Return flattened chain with burn-in discarded."""
        return self.sampler.get_chain(discard=discard, thin=thin, flat=True)

    def plot_posterior(self, discard=100, thin=10, show=True, savepath=None):
        """
        Posterior erosion rate time series: median with 68% and 95% CIs.
        """
        from matplotlib import pyplot as plt
        chain = self.get_flat_chain(discard=discard, thin=thin)
        med = np.median(chain, axis=0)
        lo1, hi1 = np.percentile(chain, [16, 84], axis=0)
        lo2, hi2 = np.percentile(chain, [2.5, 97.5], axis=0)
        # Reverse so step plots read oldest-left with invert_xaxis
        ages_plot = self.ages[::-1] / 1000
        med_r, lo1_r, hi1_r, lo2_r, hi2_r = [
            a[::-1] for a in (med, lo1, hi1, lo2, hi2)]
        fig, ax = plt.subplots()
        ax.fill_between(ages_plot[:-1], lo2_r, hi2_r,
                        step='post', alpha=0.25, color='steelblue',
                        label='95% CI')
        ax.fill_between(ages_plot[:-1], lo1_r, hi1_r,
                        step='post', alpha=0.5, color='steelblue',
                        label='68% CI')
        ax.step(ages_plot[:-1], med_r, where='post',
                color='steelblue', linewidth=2, label='Median')
        ax.set_xlabel('Age [ka]')
        ax.set_ylabel('Erosion rate [mm/yr]')
        ax.invert_xaxis()
        ax.legend()
        plt.tight_layout()
        if savepath:
            plt.savefig(savepath)
        if show:
            plt.show()
        plt.close(fig)
        return fig, ax

    def plot_predicted(self, discard=100, thin=10, n_draws=200,
                       show=True, savepath=None):
        """
        Posterior predictive: sampled model 10Be curves overlaid on data.
        """
        from matplotlib import pyplot as plt
        chain = self.get_flat_chain(discard=discard, thin=thin)
        rng = np.random.default_rng()
        idx = rng.choice(len(chain),
                         size=min(n_draws, len(chain)), replace=False)
        # Pre-compute all curves before opening a figure window
        curves = []
        for rates in chain[idx]:
            self._run_forward(rates)
            curves.append((
                self.ce.model_io['Age [yr BP]'].values.copy(),
                self.ce.model_io[
                    'Modeled surface [10Be] [atoms/g]'].values.copy()
            ))
        fig, ax = plt.subplots()
        for age_vals, conc_vals in curves:
            ax.plot(age_vals / 1000, conc_vals,
                    color='steelblue', alpha=0.05, linewidth=0.5)
        ax.errorbar(self.ce.crn_data['Age [yr BP]'] / 1000,
                    self.ce.crn_data['10Be Concentration [atoms/g]'],
                    yerr=self.ce.crn_data['10Be SD [atoms/g]'],
                    fmt='o', color='k', elinewidth=1.5, capsize=3,
                    label='Data')
        ax.set_xlabel('Age [ka]')
        ax.set_ylabel('[¹⁰Be] [atoms/g]')
        ax.invert_xaxis()
        ax.legend()
        plt.tight_layout()
        if savepath:
            plt.savefig(savepath)
        if show:
            plt.show()
        plt.close(fig)
        return fig, ax

    def save_chain(self, filepath, discard=100, thin=10):
        """
        Save the flattened posterior chain to a CSV file.
        Columns are named by erosion rate interval (e.g. eps_100.0ka-80.0ka).
        """
        chain = self.get_flat_chain(discard=discard, thin=thin)
        col_names = [
            'eps_{:.1f}ka-{:.1f}ka_mm_yr'.format(
                self.ages[k] / 1000, self.ages[k + 1] / 1000)
            for k in range(self.n_params)
        ]
        pd.DataFrame(chain, columns=col_names).to_csv(filepath, index=False)

    def save_summary_stats(self, filepath, discard=100, thin=10):
        """
        Save per-interval posterior summary (median, 68%, 95% CIs) to CSV.
        """
        chain = self.get_flat_chain(discard=discard, thin=thin)
        rows = []
        for k in range(self.n_params):
            col = chain[:, k]
            rows.append({
                'age_old_yr_BP': int(self.ages[k]),
                'age_young_yr_BP': int(self.ages[k + 1]),
                'median_mm_yr': np.median(col),
                'p16_mm_yr': np.percentile(col, 16),
                'p84_mm_yr': np.percentile(col, 84),
                'p2p5_mm_yr': np.percentile(col, 2.5),
                'p97p5_mm_yr': np.percentile(col, 97.5),
            })
        pd.DataFrame(rows).to_csv(filepath, index=False)

    def plot_summary(self, discard=100, thin=10, n_draws=200,
                     show=True, savepath=None):
        """
        Two-panel summary figure with a shared age axis.
        (a) Posterior predictive [10Be] curves overlaid on data.
        (b) Posterior erosion rate time series with 68% and 95% CIs.
        """
        from matplotlib import pyplot as plt

        chain = self.get_flat_chain(discard=discard, thin=thin)
        rng = np.random.default_rng()
        idx = rng.choice(len(chain),
                         size=min(n_draws, len(chain)), replace=False)

        # Pre-compute model curves before opening any figure window
        curves = []
        for rates in chain[idx]:
            self._run_forward(rates)
            curves.append((
                self.ce.model_io['Age [yr BP]'].values.copy(),
                self.ce.model_io[
                    'Modeled surface [10Be] [atoms/g]'].values.copy()
            ))

        # Posterior statistics for erosion rate panel
        med = np.median(chain, axis=0)
        lo1, hi1 = np.percentile(chain, [16, 84], axis=0)
        lo2, hi2 = np.percentile(chain, [2.5, 97.5], axis=0)
        ages_plot = self.ages[::-1] / 1000
        med_r, lo1_r, hi1_r, lo2_r, hi2_r = [
            a[::-1] for a in (med, lo1, hi1, lo2, hi2)]

        fig, (ax_conc, ax_eros) = plt.subplots(
            2, 1, figsize=(9, 7), sharex=True, constrained_layout=True)

        # --- Panel (a): posterior predictive [10Be] ---
        for age_vals, conc_vals in curves:
            ax_conc.plot(age_vals / 1000, conc_vals,
                         color='steelblue', alpha=0.05, linewidth=0.5)
        ax_conc.errorbar(
            self.ce.crn_data['Age [yr BP]'] / 1000,
            self.ce.crn_data['10Be Concentration [atoms/g]'],
            yerr=self.ce.crn_data['10Be SD [atoms/g]'],
            fmt='o', color='k', elinewidth=1.5, capsize=3,
            zorder=5, label='Data')
        ax_conc.set_ylabel('[¹⁰Be] [atoms/g]')
        ax_conc.legend(loc='upper right', fontsize=9)
        ax_conc.text(0.02, 0.95, '(a)', transform=ax_conc.transAxes,
                     va='top', fontsize=11, fontweight='bold')

        # --- Panel (b): posterior erosion rates ---
        ax_eros.fill_between(ages_plot[:-1], lo2_r, hi2_r,
                             step='post', alpha=0.25, color='steelblue',
                             label='95% CI')
        ax_eros.fill_between(ages_plot[:-1], lo1_r, hi1_r,
                             step='post', alpha=0.5, color='steelblue',
                             label='68% CI')
        ax_eros.step(ages_plot[:-1], med_r, where='post',
                     color='steelblue', linewidth=2, label='Median')
        ax_eros.set_xlabel('Age [ka]')
        ax_eros.set_ylabel('Erosion rate [mm/yr]')
        ax_eros.legend(loc='upper right', fontsize=9)
        ax_eros.text(0.02, 0.95, '(b)', transform=ax_eros.transAxes,
                     va='top', fontsize=11, fontweight='bold')

        ax_eros.invert_xaxis()  # shared axis: inverts both panels
        if savepath:
            plt.savefig(savepath, dpi=150, bbox_inches='tight')
        if show:
            plt.show()
        plt.close(fig)
        return fig, (ax_conc, ax_eros)
    

####################################
# ANALYTICAL SEQUENTIAL INVERSION  #
####################################

class CosmicAnalytical(object):
    """
    Sequential analytical inversion for cosmogenic 10Be erosion rate reconstruction.

    Uses data ages as the time grid. The oldest datum yields a steady-state
    erosion rate (ε = P₀Λ/C); each subsequent datum yields a rate via a 1-D
    root-find on the exact ODE solution:

        C₂ = C₁ · exp(−ε·Δt/Λ) + (P₀Λ/ε) · [1 − exp(−ε·Δt/Λ)]

    where ε appears implicitly (in both the exponential and the P₀Λ/ε
    prefactor), so brentq is required.

    Nearby samples are collapsed into inverse-variance-weighted clusters before
    solving (the sequential method assumes one datum per time step).

    Uncertainty is propagated by Monte Carlo: concentrations are perturbed by
    their 1-sigma errors and the full sequential solve is repeated for each
    replicate.

    Requires: scipy
    """

    def __init__(self, P0, attenuation_length, crn_data,
                 erosion_rate_min=0.01, erosion_rate_max=5.0, cluster_dt=500):
        """
        :param P0: Surface 10Be production rate [atoms/g/yr].
        :param attenuation_length: Attenuation length [m].
        :param crn_data: Path to CSV or DataFrame with columns
                         'Age [yr BP]', '10Be Concentration [atoms/g]',
                         '10Be SD [atoms/g]'.
        :param erosion_rate_min: Lower bound for root-find [mm/yr].
        :param erosion_rate_max: Upper bound for root-find [mm/yr].
        :param cluster_dt: Samples within this many years of each other are
                           merged into an inverse-variance-weighted mean before
                           solving. Default 500 yr.
        """
        self.P0 = P0
        self.Lambda = attenuation_length
        self.erosion_rate_min = erosion_rate_min
        self.erosion_rate_max = erosion_rate_max
        self.cluster_dt = cluster_dt

        if isinstance(crn_data, str):
            self.crn_data = pd.read_csv(crn_data)
        else:
            self.crn_data = crn_data.copy()

        self.crn_data = (self.crn_data
                         .sort_values('Age [yr BP]', ascending=False)
                         .reset_index(drop=True))

        (self.ages_clustered,
         self.conc_clustered,
         self.sigma_clustered) = self._cluster()

    def _cluster(self):
        """
        Collapse samples within cluster_dt [yr] of each other into
        inverse-variance-weighted means. Greedy forward grouping: the oldest
        sample in each group is the reference point.
        """
        ages   = self.crn_data['Age [yr BP]'].values.astype(float)
        concs  = self.crn_data['10Be Concentration [atoms/g]'].values.astype(float)
        sigmas = self.crn_data['10Be SD [atoms/g]'].values.astype(float)

        used = np.zeros(len(ages), dtype=bool)
        ages_c, concs_c, sigmas_c = [], [], []

        for i in range(len(ages)):
            if used[i]:
                continue
            group = np.where((ages[i] - ages <= self.cluster_dt) & ~used)[0]
            used[group] = True
            w = 1.0 / sigmas[group] ** 2
            ages_c.append(np.average(ages[group], weights=w))
            concs_c.append(np.average(concs[group], weights=w))
            sigmas_c.append(1.0 / np.sqrt(w.sum()))

        return np.array(ages_c), np.array(concs_c), np.array(sigmas_c)

    def _forward_step(self, C_prev, epsilon_mm, dt):
        """
        Exact ODE solution for surface [10Be] after one interval.

        C_prev    : concentration at the older time [atoms/g]
        epsilon_mm: erosion rate during the interval [mm/yr]
        dt        : interval duration [yr], positive

        Returns modeled [10Be] at the younger boundary.
        """
        eps_m = epsilon_mm / 1e3                      # mm/yr → m/yr (matches CosmicErosion)
        C_ss  = self.P0 * self.Lambda / eps_m         # steady-state [atoms/g]
        decay = np.exp(-eps_m * dt / self.Lambda)
        return C_prev * decay + C_ss * (1.0 - decay)

    def _solve_one(self, C_prev, C_obs, dt):
        """
        Root-find for erosion rate [mm/yr] given C_prev, C_obs, and interval dt.

        Returns the root within [erosion_rate_min, erosion_rate_max].  If the
        bracket does not straddle zero (concentration out of achievable range),
        returns the nearest bound and issues a warning.
        """
        from scipy.optimize import brentq

        def residual(eps):
            return self._forward_step(C_prev, eps, dt) - C_obs

        f_lo = residual(self.erosion_rate_min)
        f_hi = residual(self.erosion_rate_max)

        if f_lo * f_hi > 0:
            warnings.warn(
                'CosmicAnalytical: observed concentration {:.3g} atoms/g is '
                'outside the achievable range for this interval; clamping to '
                'nearest erosion rate bound.'.format(C_obs))
            return (self.erosion_rate_min if abs(f_lo) < abs(f_hi)
                    else self.erosion_rate_max)

        return brentq(residual, self.erosion_rate_min, self.erosion_rate_max,
                      xtol=1e-8, rtol=1e-8)

    def solve(self, ages=None, concs=None):
        """
        Sequential root-find for erosion rates.

        Uses clustered data by default; pass alternate arrays for bootstrap use.

        erosion_rates[0]  steady-state rate before the oldest datum (ε = P₀Λ/C₁)
        erosion_rates[i]  rate during the interval [ages[i-1], ages[i]], i ≥ 1

        Returns (erosion_rates [mm/yr], ages [yr BP]).
        """
        point_estimate = (ages is None)
        if ages is None:
            ages  = self.ages_clustered
            concs = self.conc_clustered

        n = len(ages)
        erosion_rates = np.zeros(n)
        C_surface     = np.zeros(n)

        # Step 0: steady-state from oldest datum
        eps0_m         = self.P0 * self.Lambda / concs[0]   # m/yr
        erosion_rates[0] = np.clip(eps0_m * 1e3,            # → mm/yr
                                   self.erosion_rate_min,
                                   self.erosion_rate_max)
        C_surface[0] = concs[0]

        # Steps 1..n-1: root-find, older → younger
        for i in range(1, n):
            dt = ages[i - 1] - ages[i]                      # positive [yr]
            erosion_rates[i] = self._solve_one(C_surface[i - 1], concs[i], dt)
            C_surface[i]     = self._forward_step(C_surface[i - 1],
                                                  erosion_rates[i], dt)

        # Only store point-estimate results; bootstrap calls must not overwrite.
        if point_estimate:
            self._last_erosion_rates = erosion_rates
            self._last_C_surface     = C_surface
        return erosion_rates, ages

    def bootstrap(self, n_boot=5000):
        """
        Monte Carlo uncertainty propagation.

        Perturbs clustered concentrations by their 1-sigma errors (Gaussian),
        repeats the sequential solve, and stores the ensemble as
        self.boot_rates (shape: n_boot × n_clustered).
        """
        rng = np.random.default_rng()
        n   = len(self.ages_clustered)
        boot_rates = np.zeros((n_boot, n))

        for k in range(n_boot):
            concs_k = (self.conc_clustered
                       + rng.standard_normal(n) * self.sigma_clustered)
            concs_k = np.clip(concs_k, 1.0, None)
            rates_k, _ = self.solve(ages=self.ages_clustered, concs=concs_k)
            boot_rates[k] = rates_k

        self.boot_rates = boot_rates
        return boot_rates

    def summary_stats(self):
        """
        DataFrame of median and 16th/84th percentile erosion rates per cluster.
        Requires bootstrap() to have been called first.
        """
        if not hasattr(self, 'boot_rates'):
            raise RuntimeError('Call bootstrap() before summary_stats().')
        median = np.median(self.boot_rates, axis=0)
        lo1, hi1 = np.percentile(self.boot_rates, [16, 84], axis=0)
        lo2, hi2 = np.percentile(self.boot_rates, [2.5, 97.5], axis=0)
        return pd.DataFrame({
            'age_yr_BP':    self.ages_clustered,
            'median_mm_yr': median,
            'p16_mm_yr':    lo1,
            'p84_mm_yr':    hi1,
            'p2p5_mm_yr':   lo2,
            'p97p5_mm_yr':  hi2,
        })

    def save_summary(self, filepath):
        """Save summary statistics to CSV."""
        self.summary_stats().to_csv(filepath, index=False)

    def plot_summary(self, show=True, savepath=None):
        """
        Two-panel summary figure.
        (a) Modeled [10Be] trajectory overlaid on observed data.
        (b) Erosion rate time series with 1-sigma bootstrap uncertainty band.

        Requires solve() and bootstrap() to have been called first.
        """
        from matplotlib import pyplot as plt

        if not hasattr(self, '_last_C_surface'):
            self.solve()

        stats  = self.summary_stats()
        ages_ka = stats['age_yr_BP'].values / 1000
        median  = stats['median_mm_yr'].values
        lo1     = stats['p16_mm_yr'].values
        hi1     = stats['p84_mm_yr'].values
        lo2     = stats['p2p5_mm_yr'].values
        hi2     = stats['p97p5_mm_yr'].values

        fig, (ax_conc, ax_eros) = plt.subplots(
            2, 1, figsize=(9, 7), sharex=True, constrained_layout=True)

        # --- Panel (a): modeled vs observed [10Be] ---
        ax_conc.plot(self.ages_clustered / 1000, self._last_C_surface,
                     color='steelblue', linewidth=2, label='Model (clustered)')
        ax_conc.errorbar(
            self.crn_data['Age [yr BP]'] / 1000,
            self.crn_data['10Be Concentration [atoms/g]'],
            yerr=self.crn_data['10Be SD [atoms/g]'],
            fmt='o', color='k', elinewidth=1.5, capsize=3,
            zorder=5, label='Data')
        ax_conc.set_ylabel('[¹⁰Be] [atoms/g]')
        ax_conc.legend(loc='upper right', fontsize=9)
        ax_conc.text(0.02, 0.95, '(a)', transform=ax_conc.transAxes,
                     va='top', fontsize=11, fontweight='bold')

        # --- Panel (b): erosion rate with bootstrap uncertainty ---
        # rate[0] is steady-state before ages[0]; rate[i] spans ages[i-1]→ages[i].
        # Build step arrays for rates[1:] (one horizontal segment per interval).
        n = len(ages_ka)
        x_step = np.empty(2 * (n - 1))
        x_step[::2]  = ages_ka[:-1]   # older boundary of each interval
        x_step[1::2] = ages_ka[1:]    # younger boundary

        def make_y_step(vals):
            y = np.empty(2 * (n - 1))
            y[::2]  = vals[1:]
            y[1::2] = vals[1:]
            return y

        ax_eros.fill_between(x_step, make_y_step(lo2), make_y_step(hi2),
                             alpha=0.25, color='steelblue', label='95% CI')
        ax_eros.fill_between(x_step, make_y_step(lo1), make_y_step(hi1),
                             alpha=0.5, color='steelblue', label='68% CI')
        ax_eros.plot(x_step, make_y_step(median),
                     color='steelblue', linewidth=2, label='Median')
        # Steady-state rate (before oldest datum) shown as a marker
        ax_eros.scatter([ages_ka[0]], [median[0]], color='steelblue',
                        zorder=5, label='Steady-state (oldest)')
        ax_eros.set_xlabel('Age [ka BP]')
        ax_eros.set_ylabel('Erosion rate [mm/yr]')
        ax_eros.legend(loc='upper right', fontsize=9)
        ax_eros.text(0.02, 0.95, '(b)', transform=ax_eros.transAxes,
                     va='top', fontsize=11, fontweight='bold')

        ax_eros.invert_xaxis()

        if savepath:
            fig.savefig(savepath, dpi=150, bbox_inches='tight')
        if show:
            plt.show()
        plt.close(fig)
        return fig, (ax_conc, ax_eros)

    def plot_comparison(self, color_ss='tab:red', show=True, savepath=None):
        """
        Compare the sequential analytical erosion rates against the naive
        steady-state estimate (ε_ss = P₀Λ/C) at each raw data point.

        Steady-state rates are shown as semi-transparent error bars.
        Steady-state uncertainty: δε_ss = ε_ss · (σ_C / C)  (first-order propagation).

        Requires solve() and bootstrap() to have been called first.
        """
        from matplotlib import pyplot as plt

        if not hasattr(self, '_last_C_surface'):
            self.solve()

        stats   = self.summary_stats()
        ages_ka = stats['age_yr_BP'].values / 1000
        median  = stats['median_mm_yr'].values
        lo1     = stats['p16_mm_yr'].values
        hi1     = stats['p84_mm_yr'].values
        lo2     = stats['p2p5_mm_yr'].values
        hi2     = stats['p97p5_mm_yr'].values

        # Steady-state estimate at every raw data point
        raw_ages   = self.crn_data['Age [yr BP]'].values / 1000
        raw_conc   = self.crn_data['10Be Concentration [atoms/g]'].values
        raw_sigma  = self.crn_data['10Be SD [atoms/g]'].values
        eps_m      = self.P0 * self.Lambda / raw_conc          # m/yr
        eps_ss     = eps_m * 1e3                               # mm/yr
        eps_ss     = np.clip(eps_ss, self.erosion_rate_min,
                             self.erosion_rate_max)
        eps_ss_err = eps_ss * (raw_sigma / raw_conc)           # 1σ, first-order

        # Step function for analytical rates
        n = len(ages_ka)
        x_step = np.empty(2 * (n - 1))
        x_step[::2]  = ages_ka[:-1]
        x_step[1::2] = ages_ka[1:]

        def make_y_step(vals):
            y = np.empty(2 * (n - 1))
            y[::2]  = vals[1:]
            y[1::2] = vals[1:]
            return y

        fig, ax = plt.subplots(figsize=(9, 4), constrained_layout=True)

        # Steady-state: error bars at each raw data point
        ax.errorbar(raw_ages, eps_ss, yerr=eps_ss_err,
                    fmt='o', color=color_ss, alpha=0.6,
                    elinewidth=1.5, capsize=3, zorder=3,
                    label='Steady-state (per datum)')

        # Analytical: CIs and median step function
        ax.fill_between(x_step, make_y_step(lo2), make_y_step(hi2),
                        alpha=0.25, color='steelblue', label='Analytical 95% CI')
        ax.fill_between(x_step, make_y_step(lo1), make_y_step(hi1),
                        alpha=0.5,  color='steelblue', label='Analytical 68% CI')
        ax.plot(x_step, make_y_step(median),
                color='steelblue', linewidth=2, label='Analytical median')
        ax.scatter([ages_ka[0]], [median[0]], color='steelblue',
                   zorder=5)

        ax.set_xlabel('Age [ka BP]')
        ax.set_ylabel('Erosion rate [mm/yr]')
        ax.invert_xaxis()
        ax.legend(loc='upper right', fontsize=9)

        if savepath:
            fig.savefig(savepath, dpi=150, bbox_inches='tight')
        if show:
            plt.show()
        plt.close(fig)
        return fig, ax


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
            
