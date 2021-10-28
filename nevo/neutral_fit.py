import os
import scipy
import numpy as np
import pandas as pd
from datetime import datetime
from lmfit import Parameters, Model, fit_report 
from scipy.stats import beta 
from statsmodels.stats.proportion import proportion_confint 
from nevo.neutral_fit_utils import beta_cdf, subsample
from nevo.neuplot import neufit_plot

def nevo_pipeline(output_filename, dataset_type, custom_filename, 
                  norm_graph = True, colored_graph = True, non_neutral = True, 
                  non_save = False, full_non_neutral = False):
    
    '''Calls all functions needed to create neutral model 
    
    Written by: Caitlin Guccione, 08-25-2021
    
    Parameters
    ----------
    output_filename: str
        Name/nickname of dataset (ex. 'combined'). Will be incorperated into 
        filename of all nevo outputs. 
    dataset_type: str
        Type of data importing. Currently just ('hutchKraken', 'gregTCGA'), 
        in ToDo, make more specific to be used more generally.
    custom_filename: str
        An extra random varible depedatn on dataset
            - hutchKraken : the name of the biom file ex.'combined_biome'
    norm_graph: bool, optional
        If 'False', does NOT print or save the neutral evolution graph without 
        any custom bacteria coloring
    colored_graph: bool, optional
        If 'False', does NOT prints or save the neutral evolution graph with 
        custom bacteria coloring
    non_neutral: bool, optional
        If 'False', does NOT print or save the most non-neutral microbes into 
        a .csv file
    non_save: bool, optional
        If 'True', will only prints and does NOT save: 
        the neutral evolution graph without any custom bacteria coloring, 
        the neutral evolution graph with custom bacteria coloring, and 
        the most non-neutral microbes csv file.
        This was intened for testing purposes. 
    full_non_neutral: bool, optional
        If 'True', will create an additonal csv file about which points are 
        neutral, and how far off the curve each point is. The file csv file
        will contain: otu_id, mean_abundance, occurance, Kingdom, Phylum,
        Class, Order, Family, Genus, Species, predicted_occurence, 
        lower_conf_int, and upper_conf_int. This data can be used for custom
        coloring of the neutral evolution graph. 

    TODO
    ----
    - Change dataset_type from 'hutchKraken' / 'gregTCGA' to biom input files, 
    biom input files missing taxonomy, ...
    - Remove the 'custom_filename' varible and create another funciton that 
    handles the current issue 
    - Remove everything from '#Convert data from biom to csv files for Neufit'
    and below because it is too specific
    - Determine how to not create Neufit text file that ends up getting 
    deleted
    
    '''
    
    #Convert data from biom to csv files for Neufit
    if dataset_type == 'hutchKraken':
        fnData, fnTaxonomy = biom2data_tax(hutchKrakenAlex_biom, 
                                           custom_filename, output_filename)
    elif dataset_type == 'TCGA_WGS':
        #Determine which biom file for the run
        if output_filename == 'normal':
            biomFilename ='116640_feature-table-TCGA-WGS-STN-ESCA-HNSC.biom'
        elif output_filename == 'cancer':
            biomFilename ='116639_feature-table-TCGA-WGS-PT-ESCA-HNSC.biom'
        
        #Turn biome file into _data and _tax files
        fnD_hn, fnD_e, fnT_hn, fnT_e, fnD_e_scc, fnT_e_scc, fnD_e_eac, 
        fnT_e_eac = biom_addMetaTax_customTCGAehn(tcgaEhnWGSgreg_, 
                                                  biomFilename, 
                                                  output_filename)
        
        #Choose the correct _data file, esoph or head and neck
        if custom_filename == 'e':
            fnData = fnD_e
            fnTaxonomy = fnT_e
            output_filename = 'esophagus_' + output_filename
        elif custom_filename == 'hn':
            fnData = fnD_hn
            fnTaxonomy = fnT_hn
            output_filename = 'headNeck_' + output_filename
        elif custom_filename == 'e_scc':
            fnData = fnD_e_scc
            fnTaxonomy = fnT_e_scc
            output_filename = 'esophagus_squamousCellCarcinoma_' + \ 
            output_filename
        elif custom_filename == 'e_eac':
            fnData = fnD_e_eac
            fnTaxonomy = fnT_e_eac
            output_filename = 'esophagus_adenocarcinoma_' + \ 
            output_filename
        
        
    #Run Neufit
    occurr_freqs, n_reads, n_samples, r_square, beta_fit, 
    file_header = main_neufit(output_filename, dataset_type, fnData, 
                              fnTaxonomy, full_non_neutral)
    
    #Neufit Plotting and Non-neutral Outline
    if norm_graph == True: #Neutral evolution graph, no color
        nc_fn = neufit_plot(occurr_freqs, n_reads, n_samples, r_square, 
                            beta_fit, file_header)
        if non_save == False:
            save_plot(nc_fn)
    if colored_graph == True:#Neutral evolution graph with colors
        cc_fn = custom_color_plot(occurr_freqs, n_reads, n_samples, 
                                  r_square, beta_fit, file_header)
        if non_save == False:
            save_plot(cc_fn)
    if non_neutral == True:
        non_neutral_outliers(file_header, occurr_freqs, dataset_type, 
                             non_save)
        
    #Easier to delete the Neufit text file then to not create it
    if non_save == True:
        neufit_fn= str(file_header) + ".txt"
        os.remove(neufit_fn)


def neufit(output_filename, dataset_type, _data_filename, _taxonomy_filename, 
           full_non_neutral = False, arg_ignore_level = 0, 
           arg_rarefaction_level = 0):
    
    '''Fits a neutral community model to species abundances
    
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
    output_filename: str
        Name/nickname of dataset (ex. 'combined'). Will be incorperated into 
        filename of all nevo outputs.     
    dataset_type: str
        Type of data importing. Currently just ('hutchKraken', 'gregTCGA'), 
        in ToDo, make more genearl
    _data_filename: str, path
        The path to []_data.csv file; often an OTU abudance table. 
    _taxonomy_filename: str, path
        The path to []_taxonomy.csv; corresponding taxonomic information.
    arg_ignore_level: int, optional
        Ignores OTUs below this abudance threshold; default is to use all 
        OTUs regardless of abudance threshold. Value must be non-negative.
    arg_rarefaction_level: int, optional
        Sets the rarefaction level. Leaving the default of 0 changes 
        this value to the highest possible uniform read depth. 
    
    Returns
    -------
    file_header: str, path
        Filepath for all nevo outputs. Includes path, data nickname and 
        time stamp.
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

    Notes
    -----
    The orginal code for this function comes from the neufit.py file
    (https://github.com/misieber/neufit/blob/master/neufit.py) and credit
    is given to Michael Sieber for this function. The following minor
    changes were made in the function below in order to make the 
    code functional in the pipeline:
        - Converted code in neufit.py into a function by changing all argv 
        inputs into function inputs.
            - args.data_file.name -> data_file_name
            - args.data_file <open file 'data_.csv', mode 'r' 
            at 0x120b28c00> -> data_file = open(data_file_name, "r")
            - args.ignore_level 0 -> arg_ignore_level = 0
            - args.rarefaction_level 0 -> arg_rarefaction_level = 0
            - args.taxonomy_file <open file '_taxonomy.csv', mode 'r' 
            at 0x120b28c90> -> tax_file = open(tax_file_name, "r")
        - Changed was all syntax from Python2 to Python3.
        - Changed all print statments to write to file statements
            - Made any necessary edits to allow data to be printed to 
            file instead of terminal
        - Output pandas df that is created to be used in future steps
    
    TODO
    ----
    - Change dataset_type from 'hutchKraken' / 'gregTCGA' to biom 
    input files, biom input files missing taxonomy, ...
    - Remove the 'custom_filename' varible and create another funciton that 
    handles the current issue 
    - Add test cases for arg_ignore_level and arg_arefaction_level to insure
    non-negative input 
    
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
    
    ##Added by Caitlin ~ Push output to file instead of printing to screen
    
    #Grab and format data/time
    time = datetime.time(datetime.now())
    date = datetime.date(datetime.now())
    h,s = str(time).split(".") #Split string into hours/min and sec
    
    #Create file_header which holds the path for all future NEvo outpus
    file_header = str(neufit_output_path) + "/" + str(dataset_type) + \
        '/' + str(output_filename) + '/' + str(output_filename) + '_' + \
        str(date) + "_" + str(h)
    
    #Creates directory for all NEvo Outputs if it doesn't exist 
    dir_name = str(neufit_output_path) + "/" + str(dataset_type) + \
        '/' + str(output_filename)
    os.makedirs(dir_name, exist_ok=True)
    
    #Create and open file for Neufit Output txt file
    fn= str(file_header) + ".txt"
    file = open(fn, 'w')
    
    #Print statments with important info
    print("Running dataset: " + str(dataset_type) + "Category:" + \
          str(output_filename) + '\n')
    ##
    
    # Writes dataset info output file, calculates and writes the 
    # number of samples/ reads in the file
    file.write('Corresponding csv file: ' + _data_filename + '\n')
    abundances = pd.read_table(_data_filename, header=0, 
                               index_col=0, sep='\t').astype(int)
    abundances = abundances[abundances.sum(1) > arg_ignore_level]
    file.write ('Dataset contains ' + str(abundances.shape[1]) + \
                ' samples (sample_id, reads): \n')
    
    ##Caitlin
    # The following loop is used instead of 'print abundances.sum(0)' 
    # so that it can be written to a file
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
    file.write ('(' + str(arg_rarefaction_level) + \
                ' reads per sample) \n')

    # Optionally subsample the abundance table, unless all samples 
    # already have the required uniform read depth
    if not all(n_reads == arg_rarefaction_level for n_reads in abundances.sum(0)):
        abundances = subsample(abundances, arg_rarefaction_level)
        abundances = abundances[abundances.sum(1) > 0]

    # Dataset shape
    n_otus, n_samples = abundances.shape
    n_reads = arg_rarefaction_level

    file.write ('fitting neutral expectation to dataset with ' + \
                str(n_samples) + ' samples and ' + str(n_otus) + \
                ' otus \n \n')
    # Calculate mean relative abundances and occurrence frequencies
    mean_relative_abundance = (1.0*abundances.sum(1))/n_reads/n_samples
    occurrence_frequency = (1.0*np.count_nonzero(abundances, 
                                                 axis=1))/n_samples
    
    occurr_freqs = pd.DataFrame(mean_relative_abundance, 
                                columns=['mean_abundance'])
    if dataset_type == 'TCGA_WGS':#This changes name of first column
        occurr_freqs.index.name = 'gOTU' 
    else:
        occurr_freqs.index.name = 'otu_id'
    occurr_freqs['occurrence'] = occurrence_frequency
    occurr_freqs = occurr_freqs.sort_values(by=['mean_abundance'])

    # Join with taxonomic information (optional)
    if _taxonomy_filename != None: #Changed <> to !=
        if dataset_type == 'TCGA_WGS':
             taxonomy = pd.read_table(_taxonomy_filename, header=0, 
                                      index_col=1, sep='\t')
        else:
            taxonomy = pd.read_table(_taxonomy_filename, header=0, 
                                     index_col=0, sep='\t')
        occurr_freqs = occurr_freqs.join(taxonomy)
        
    # Fit the neutral model
    params = Parameters()
    params.add('N', value=n_reads, vary=False)
    params.add('m', value=0.5, min=0.0, max=1.0)
    beta_model = Model(beta_cdf)
    beta_fit = beta_model.fit(occurr_freqs['occurrence'], params, 
                              p=occurr_freqs['mean_abundance'])

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
        pd.concat((above, below)).to_csv(str(file_header) + \
                                         '_FullNonNeutral.csv')
    
    file.close()
    
    return(occurr_freqs, n_reads, n_samples, r_square, beta_fit, file_header)