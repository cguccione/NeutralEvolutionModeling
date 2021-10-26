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

import os
import scipy
import numpy as np
import pandas as pd
from datetime import datetime
from lmfit import Parameters, Model, fit_report 
from scipy.stats import beta 
from statsmodels.stats.proportion import proportion_confint 
from neuevo.neufit_utils import beta_cdf, subsample
from neuevo.neuplot import neufit_plot

'''
All the code in the function below was written by the orginal Neufit authors, not me and is from the neufit.py file.
neufit.py: https://github.com/misieber/neufit/blob/master/neufit.py

Here are the following changes I made:
- Turned it into a function instead of using argv as input
    -There are specific notes on what I changed in the 'Notes: Neufit Modifyed Args Input/Output Details' section.
-Changed was all syntax from Python2 to Python3.
-Changed all print statments to write to file statements
    - Made any necissary edits to allow data to be printed to file instead of terminal 
'''


def neufit(fnData, fnTaxonomy, output_filename,
           dataset_type, custom_filename,
           norm_graph, colored_graph, non_neutral,
           non_save = False, full_non_neutral = False):
    '''
    Inputs:    output_filename : the name which will be at the front of all the files; ex. 'combined'
               dataset_type : the dataset type we have from here: ('hutchKraken', 'gregTCGA')
               custom_filename : depedant on the dataset:
                       - hutchKraken : the name of the biom file ex.'combined_biome'
               
               norm_graph : True/False : Prints and saves the neutral evolution graph without any coloring
               colored_graph : True/False : Prints and saves the neutral evolution graph without any coloring
               non_neutral : True/False : Prints and saves the most non-neutral microbes in csv file               
               
               Default inputs to be changed mainly for testing purposes:
               non_save : True/False, False = default, the following will not be SAVED just printed to the screen
                           - * This is intened for testing only! *
                           - norm_graph, colored_graph, non_neutral_csv 
                           - **Took away this feature ** [It will still create [name]_data.csv and [name]_taxonomy.csv but will delete them after running]
               full_non_neutral: True/False : Creates a csv file from orginal Neufit program with information about what is neutral and how far off the curve each point is
                   - csv will have: otu_id, mean_abundance, occurrence, Kingdom, Phylum, Class, Order, Family, Genus, Species, predicted_occurrence, lower_conf_int, upper_conf_int

              *The main function for the pipeline which calls all other functions*
    
    '''

    #Run Neufit
    occurr_freqs, n_reads, n_samples, r_square, beta_fit, file_header = main_neufit(output_filename, 
                                                                                    dataset_type,
                                                                                    fnData, fnTaxonomy,
                                                                                    full_non_neutral)
    
    #Neufit Plotting and Non-neutral Outline
    if norm_graph == True: #Neutral evolution graph, no color
        nc_fn = neufit_plot(occurr_freqs, n_reads, n_samples, r_square, beta_fit, file_header)
        if non_save == False:
            save_plot(nc_fn)
    if colored_graph == True:#Neutral evolution graph with colors
        cc_fn = custom_color_plot(occurr_freqs, n_reads, n_samples, r_square, beta_fit, file_header)
        if non_save == False:
            save_plot(cc_fn)
    if non_neutral == True:
        non_neutral_outliers(file_header, occurr_freqs, dataset_type, non_save)
        
    #Easier to delete the Neufit text file then to not create it
    if non_save == True:
        neufit_fn= str(file_header) + ".txt"
        os.remove(neufit_fn)
    
    #Optional cleanup step with non_save to remove taxonomy and data files - these just overwrite themseleves so not usally a big issue
    '''
    if non_save == True: #Delete [name]_data.csv and [name]_taxonomy.csv to reduce cluter
        os.remove(fnData)
        os.remove(fnTaxonomy)
    '''


def main_neufit(output_filename, dataset_type, _data_filename, _taxonomy_filename, full_non_neutral = False, arg_ignore_level = 0, arg_rarefaction_level = 0):
    '''Inputs:  output_filename : the name which will be at the front of all the files; ex. 'combined'
                dataset_type : the dataset type we have from here: ('hutchKraken', 'gregTCGA')
                _data_filename = path of []_data.csv file needed for Neufit to run
                _taxonomy_filename = path of []_taxonomy.csv file needed for Neufit to run
                arg_ignore_level = 0 ; default set from orginal Neufit program
                arg_rarefaction_level = 0; default set from orginal Neufit program 
       Outputs: file_header = file path for neufit outpus: The path+dataGroup name + date stamp to be used for all parts of Neufit run: ex. /home/cguccion/NeutralEvolutionModeling/neufit_output/hutchKrakenAlex_combined_2021-08-26_13:24:22  
                occurr_freqs = pandas df that Neufit created in the orginal program but didn't specifcally output orginally
                    - Headers of csv: otu_id, mean_abundance, occurrence, Kingdom, Phylum, Class, Order, Family, Genus, Species, predicted_occurrence, lower_conf_int, upper_conf_int
                    - I used this occur_freqs df in order to figure out which species is the most non-neutral
                    - I believe this is what Neufit uses to physical plot the dots on the graph as well 
                n_reads = Number of reads (from orginal program)
                n_samples = Number of smaples (from orginal program)
                r_square = R^2 value (from orginal program)
                beta_fit = stats on the preformance of the model (from orginal program)
               *Runs the main section of neufit'''
    
    ##Added by Caitlin ~ Push output to file instead of printing to screen
    
    #Grab and format data/time
    time = datetime.time(datetime.now())
    date = datetime.date(datetime.now())
    h,s = str(time).split(".") #Split  the string into  hours/min and seconds
    
    #Create file_header which holds the path / location for all future Neufit outpus
    file_header = str(neufit_output_path) + "/" + str(dataset_type) + '/' + str(output_filename) + '/' + str(output_filename) + '_' + str(date) + "_" + str(h)
    
    #Creates directory for all Neufit Outputs if it doesn't already exist 
    dir_name = str(neufit_output_path) + "/" + str(dataset_type) + '/' + str(output_filename)
    os.makedirs(dir_name, exist_ok=True)
    
    #Create and open file for Neufit Output txt file
    fn= str(file_header) + ".txt"
    file = open(fn, 'w')
    
    #Print statments with important info
    print("Running dataset: " + str(dataset_type) + "Category:" + str(output_filename) + '\n')
    ##
    
    # Writes dataset info to Neufit output file + calculates and writes the number of samples/ reads in the file
    file.write('Corresponding csv file: ' + _data_filename + '\n')
    abundances = pd.read_table(_data_filename, header=0, index_col=0, sep='\t').astype(int)
    abundances = abundances[abundances.sum(1) > arg_ignore_level]
    file.write ('Dataset contains ' + str(abundances.shape[1]) + ' samples (sample_id, reads): \n')
    ##Caitlin
    #The following loop is used instead of 'print abundances.sum(0)' so that it can be written to a file
    for index, col in abundances.iteritems():
        col_sum = 0
        for i in col:
            col_sum += i
        file.write (index + '\t' + str(col_sum) + '\n')
    file.write ('\n')
    ##

    # Determine uniform read depth
    if arg_rarefaction_level == 0 or arg_rarefaction_level > max(abundances.sum(0)):
        arg_rarefaction_level = min(abundances.sum(0))
        file.write ('rarefying to highest possible uniform read depth'),
    else:
        file.write ('rarefying to custom rarefaction level'),
    file.write ('(' + str(arg_rarefaction_level) + ' reads per sample) \n')

    # Optionally subsample the abundance table, unless all samples already have the required uniform read depth
    if not all(n_reads == arg_rarefaction_level for n_reads in abundances.sum(0)):
        abundances = subsample(abundances, arg_rarefaction_level)
        abundances = abundances[abundances.sum(1) > 0]

    # Dataset shape
    n_otus, n_samples = abundances.shape
    n_reads = arg_rarefaction_level

    file.write ('fitting neutral expectation to dataset with ' + str(n_samples) + ' samples and ' + str(n_otus) + ' otus \n \n')
    # Calculate mean relative abundances and occurrence frequencies
    mean_relative_abundance = (1.0*abundances.sum(1))/n_reads/n_samples
    occurrence_frequency = (1.0*np.count_nonzero(abundances, axis=1))/n_samples
    
    occurr_freqs = pd.DataFrame(mean_relative_abundance, columns=['mean_abundance'])
    if dataset_type == 'TCGA_WGS':
        occurr_freqs.index.name = 'gOTU' #This changes the name of the first column
    else:
        occurr_freqs.index.name = 'otu_id'
    occurr_freqs['occurrence'] = occurrence_frequency
    occurr_freqs = occurr_freqs.sort_values(by=['mean_abundance'])

    # Join with taxonomic information (optional)
    if _taxonomy_filename != None: #Changed <> to !=
        if dataset_type == 'TCGA_WGS':
             taxonomy = pd.read_table(_taxonomy_filename, header=0, index_col=1, sep='\t')
        else:
            taxonomy = pd.read_table(_taxonomy_filename, header=0, index_col=0, sep='\t')
        occurr_freqs = occurr_freqs.join(taxonomy)
        
    # Fit the neutral model
    params = Parameters()
    params.add('N', value=n_reads, vary=False)
    params.add('m', value=0.5, min=0.0, max=1.0)
    beta_model = Model(beta_cdf)
    beta_fit = beta_model.fit(occurr_freqs['occurrence'], params, p=occurr_freqs['mean_abundance'])

    # Report fit statistics
    r_square = 1.0 - np.sum(np.square(occurr_freqs['occurrence'] - beta_fit.best_fit))/np.sum(np.square(occurr_freqs['occurrence'] - np.mean(occurr_freqs['occurrence'])))
    file.write (fit_report(beta_fit))
    file.write ('\n R^2 = ' + '{:1.2f}'.format(r_square))
    print(fit_report(beta_fit))
    print('\n R^2 = ' + '{:1.2f}'.format(r_square))
    print('=========================================================')

    # Adding the neutral prediction to results
    occurr_freqs['predicted_occurrence'] = beta_fit.best_fit
    occurr_freqs['lower_conf_int'], occurr_freqs['upper_conf_int'] = proportion_confint(occurr_freqs['predicted_occurrence']*n_samples, n_samples, alpha=0.05, method='wilson')

    # Save non-neutral otus (here simply determined by lying outside the confidence intervals)
    above = occurr_freqs[occurr_freqs['occurrence'] > occurr_freqs['upper_conf_int']]
    below = occurr_freqs[occurr_freqs['occurrence'] < occurr_freqs['lower_conf_int']]
    
    #Create orginal non neutral output file from Neufit
    if full_non_neutral == True:
        pd.concat((above, below)).to_csv(str(file_header) + '_FullNonNeutral.csv')
    
    file.close()
    
    return(occurr_freqs, n_reads, n_samples, r_square, beta_fit, file_header)