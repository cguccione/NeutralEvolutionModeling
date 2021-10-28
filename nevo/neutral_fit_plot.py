import os
import numpy as np
from scipy.stats import beta 
from datetime import datetime
from statsmodels.stats.proportion import proportion_confint 
from matplotlib import pyplot
from math import log10

def neufit_plot(occurr_freqs, n_reads, n_samples, r_square, beta_fit, file_header):
    '''Creates the neutral evolution png plot 
    
    Modified by: Caitlin Guccione 08-25-2021
    
    Copyright
    ---------
    Github: https://github.com/misieber/neufit
    Theroy as described in [1]_. and [2]_.
    
    Copyright (C) 2018 Michael Sieber (sieber.ecoevo.de)
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions 
    are met:

    1. Redistributions of source code must retain the above copyright 
       notice, this list of conditions and the following disclaimer.
    2. Redistributions in binary form must reproduce the above 
       copyright notice, this list of conditions and the following 
       disclaimer in the documentation and/or other materials provided 
       with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
    "AS IS" ANDANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
    Parameters
    ----------
    occurr_freqs: pandas df
        Df header: otu_id, mean_abundance, occurrence, Kingdom, Phylum, 
        Class, Order, Family, Genus, Species, predicted_occurrence, 
        lower_conf_int, upper_conf_int
    n_reads: int
        Total number of reads.
    n_samples: int
        Total number of samples.
    r_square: float
        R^2 value of the fit of data to neutral curve.
    beta_fit: lmfit.model.ModelResult object
        Holds the stats on the preformance of the model.
    file_header: str, path
        Filepath for all nevo outputs. Includes path, data nickname and 
        time stamp.
        
    Returns
    -------
    png
        A plot showing the the input data and how well it fits the neutral
        model 
    fn: str, path
        A filepath where the plot is stored. 

    Notes
    -----
    The orginal code for this function comes from the neufit.py file
    (https://github.com/misieber/neufit/blob/master/neufit.py) and credit
    is given to Michael Sieber for this function. The following minor
    changes were made in the function below in order to make the 
    code functional in the pipeline:
        - Turned the plotting section into a reusable function
        - Added a plot clearing option to avoid double keys
    
    TODO
    ----
    - Add a better definition for returns, is there a more clear 
    way to describe the plot
    - See if you can eliminate the return of filename from the plot 
    
    References
    ----------
    .. [1] Sloan, W. T., Lunn, M., Woodcock, S., Head, I. M., Nee, S. 
    and Curtis, T. P. (2006). Quantifying the roles of immigration 
    and chance in shaping prokaryote community structure. Environmental 
    Microbiology, 8:732-740. https://doi.org/10.1111/j.1462-2920.2005.00956.x
    
    .. [2] Sieber, M., Pita, L., Weiland-Bräuer, N., Dirksen, P., Wang, J., 
    Mortzfeld, B., Franzenburg, S., Schmitz, R. A., Baines, J. F., Fraune, 
    S., Hentschel, U., Schulenburg, H., Bosch, T. C. G. and Traulsen, A. 
    (2018). The Neutral Metaorganism. bioRxiv. https://doi.org/10.1101/367243
    '''
    
    pyplot.cla() #CG: Clears previous plot - to avoid double keys

    # Prepare results plot
    pyplot.xlabel('Mean relative abundance across samples', fontsize=15)
    pyplot.xscale('log')
    x_range = np.logspace(log10(min(occurr_freqs['mean_abundance'])/10),
                          0, 1000)
    pyplot.xlim(min(x_range), max(x_range))
    pyplot.xticks(fontsize=16)
    pyplot.ylabel('Occurrence frequency in samples', fontsize=15)
    pyplot.ylim(-0.05, 1.05)
    pyplot.yticks(fontsize=16)

    # Plot data points
    pyplot.plot(occurr_freqs['mean_abundance'], occurr_freqs['occurrence'],
                'o', markersize=6, fillstyle='full', color='black')

    # Plot best fit
    pyplot.plot(x_range, beta_cdf(x_range, n_reads, beta_fit.best_values['m']),
                '-', lw=5, color='darkred')
    lower, upper = proportion_confint(beta_cdf(x_range, n_reads,
                                               beta_fit.best_values['m'])*n_samples,
                                      n_samples, alpha=0.05, method='wilson')
    pyplot.plot(x_range, lower, '--', lw=2, color='darkred')
    pyplot.plot(x_range, upper, '--', lw=2, color='darkred')
    pyplot.fill_between(x_range, lower, upper, color='lightgrey')

    pyplot.text(0.05, 0.9, '$R^2 = ' + '{:1.2f}'.format(r_square) + '$',
                fontsize=16, transform=pyplot.gca().transAxes)
    pyplot.tight_layout()
    
    fn = file_header + '_NeutralFitPlot.png'
    return(fn)