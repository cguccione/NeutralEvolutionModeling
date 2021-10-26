# neufit: Fit a neutral community model to species abundances, e.g. from an OTU table
#
# For the theory behind this see Sloan et al, Environ Microbiol 2006 8:732-740.
# To run on the example simulation data: python neufit.py sim_data.csv
# To link with the mock taxonomy use the -t sim_taxonomy.csv option
#
# Copyright (C) 2018 Michael Sieber (sieber.ecoevo.de)
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Github: https://github.com/misieber/neufit

'''
All the code in the function below was written by the orginal Neufit authors, not me and is from the neufit.py file.
neufit.py: https://github.com/misieber/neufit/blob/master/neufit.py

Here are the following changes I made:
- Turned it into a function
- Added a plot clearing option to avoid double keys
'''
import os
import numpy as np
from scipy.stats import beta 
from datetime import datetime
from statsmodels.stats.proportion import proportion_confint 
from matplotlib import pyplot
from math import log10

def neufit_plot(occurr_freqs, n_reads, n_samples, r_square, beta_fit, file_header):
    '''Inputs: occurr_freqs = pandas df that Neufit created in the orginal program but didn't specifcally output orginally
                    - Headers of csv: otu_id, mean_abundance, occurrence, Kingdom, Phylum, Class, Order, Family, Genus, Species, predicted_occurrence, lower_conf_int, upper_conf_int
                    - I used this occur_freqs df in order to figure out which species is the most non-neutral
                    - I believe this is what Neufit uses to physical plot the dots on the graph as well
               n_reads = Number of reads (from orginal program)
               n_samples = Number of smaples (from orginal program)
               r_square = R^2 value (from orginal program)
               beta_fit = stats on the preformance of the model (from orginal program)
               file_header = file path for all neufit outpus: The path+dataGroup name + date stamp to be used for all parts of Neufit run: ex. /home/cguccion/NeutralEvolutionModeling/neufit_output/hutchKrakenAlex_combined_2021-08-26_13:24:22
        Outputs: ! The orginal plot that Neufit outputs (just black dots) - will print to the screen but not save when run
                 fn = default filename for neufit plot : customName_NeutralFitPlot.png -- Will often not need this
                 *Creates the plot for Neufit '''


    pyplot.cla() #Clears previous plot - to avoid double keys

    # Prepare results plot
    pyplot.xlabel('Mean relative abundance across samples', fontsize=15)
    pyplot.xscale('log')
    x_range = np.logspace(log10(min(occurr_freqs['mean_abundance'])/10), 0, 1000)
    pyplot.xlim(min(x_range), max(x_range))
    pyplot.xticks(fontsize=16)
    pyplot.ylabel('Occurrence frequency in samples', fontsize=15)
    pyplot.ylim(-0.05, 1.05)
    pyplot.yticks(fontsize=16)

    # Plot data points
    pyplot.plot(occurr_freqs['mean_abundance'], occurr_freqs['occurrence'], 'o', markersize=6, fillstyle='full', color='black')

    # Plot best fit
    pyplot.plot(x_range, beta_cdf(x_range, n_reads, beta_fit.best_values['m']), '-', lw=5, color='darkred')
    lower, upper = proportion_confint(beta_cdf(x_range, n_reads, beta_fit.best_values['m'])*n_samples, n_samples, alpha=0.05, method='wilson')
    pyplot.plot(x_range, lower, '--', lw=2, color='darkred')
    pyplot.plot(x_range, upper, '--', lw=2, color='darkred')
    pyplot.fill_between(x_range, lower, upper, color='lightgrey')

    pyplot.text(0.05, 0.9, '$R^2 = ' + '{:1.2f}'.format(r_square) + '$', fontsize=16, transform=pyplot.gca().transAxes)
    pyplot.tight_layout()

    fn = file_header + '_NeutralFitPlot.png'

    return(fn)