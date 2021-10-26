# Caitlin Guccione, 08-25-2021

import pandas as pd
from biom import load_table 
import os
import numpy as np

def non_neutral_outliers(file_header, occurr_freqs, dataset_type, non_save, threshold = 0.5):
    '''Inputs: file_header = file path for neufit outpus: The path+dataGroup name + date stamp to be used for all parts of Neufit run: ex. /home/cguccion/NeutralEvolutionModeling/neufit_output/hutchKrakenAlex_combined_2021-08-26_13:24:22  
               occurr_freqs = pandas df that Neufit created in the orginal program but didn't specifcally output orginally
                    - Headers of csv: otu_id, mean_abundance, occurrence, Kingdom, Phylum, Class, Order, Family, Genus, Species, predicted_occurrence, lower_conf_int, upper_conf_int
                    - I used this occur_freqs df in order to figure out which species is the most non-neutral
                    - I believe this is what Neufit uses to physical plot the dots on the graph as well
               threshold = 0.5; autoset to 0.5, but determines what bacteria are considered non-neutral
        Outputs: ! [name]_NonNeutralOutliers.csv
                *Creates the most Non-neutral csv file '''
    
    #Create dataframe
    standoutMicrobes = pd.DataFrame(columns = ('Difference off Neutral Model', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))
    
    #Loop and find most non-neutral microbes based upon threshold
    row_count = 0
    for i,j in occurr_freqs.iterrows():
        diff = abs(j['occurrence'] - j['predicted_occurrence'])
        if diff > threshold:
            if dataset_type == 'TCGA_WGS':
                standoutMicrobes.loc[row_count] = [diff, j['Domain'], j['Phylum'], j['Class'], j['Order'], j['Family'], j['Genus'], j['Species']]
            else:
                standoutMicrobes.loc[row_count] = [diff, j['Kingdom'], j['Phylum'], j['Class'], j['Order'], j['Family'], j['Genus'], j['Species']]
            row_count +=1
    standoutMicrobes = standoutMicrobes.sort_values(by =['Difference off Neutral Model'], ascending=False)
    
    #Display and export non-neutral microbes as csv
    print("\nTop NonNeutral Microbes")
    display(standoutMicrobes)
    print('=========================================================\n')
    
    fn = file_header + '_NonNeutral_Outliers.csv'
    
    if non_save == False:
        standoutMicrobes.to_csv(fn, sep = '\t')

def biom2data_tax(datasetName, biomFilename, finalFilename):
    '''Inputs: datasetName = filename of dataset ex. hutchKrakenAlex_biom -- found in set paths above
               biomFilename = filename of biom file found in datasetName path above
               finalFilename = filename to be used in _data.csv file for neuFit
       Outputs: fnD = File location of [name]_data.csv
                ![name]_data.csv = Creates a data.csv file to be input into Neufit located at 'neufit_input'
                fnT = File location of [name]_taxonomy.csv
                ![name]_taxonomy.csv = Creates a taxonomy.csv file to be input into Neufit located at 'neufit_input'
       * Imports biom file -> pandas dataframe -> data.csv, taxonomy.csv as desired by neufit'''
    
    #Make filename and import data
    fullFilename = datasetName + '/' + biomFilename
    featureTable = load_table(fullFilename) #https://biom-format.org/documentation/generated/biom.load_table.html
    
    
    '''
    #Create file_header which holds the path / location for all _data and _csv outputs for this data
    ###file_header = str(neufit_output_path) + "/" + str(dataset_type) + '/' + str(output_filename) + '/' + str(output_filename) + '_' + str(date) + "_" + str(h)
    file_header = str(neufit_output_path) + "/" + str(datasetName) + '/' + str(finalFilename) + '/' + str(finalFilename) + '_' + str(date) + "_" + str(h)
    
    #Creates directory for all Neufit Outputs if it doesn't already exist 
    ###dir_name = str(neufit_output_path) + "/" + str(dataset_type) + '/' + str(output_filename)
    dir_name = str(neufit_output_path) + "/" + str(datasetName) + '/' + str(finalFilename)
    os.makedirs(dir_name, exist_ok=True)
    '''
    
    #Create _data.csv
    pandas_featureTable = pd.DataFrame(featureTable.matrix_data.toarray(), featureTable.ids('observation'), featureTable.ids())
    fnD = neufit_input_path + '/' + finalFilename + '_data.csv'
    pandas_featureTable.to_csv(fnD, sep='\t')
    
    #Create _taxonomy.csv
    pandas_TaxTable = pd.DataFrame(featureTable.metadata_to_dataframe('observation'))
    pandas_TaxTable.set_axis(['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'], axis=1)
    fnT = neufit_input_path + '/' + finalFilename + '_taxonomy.csv'
    pandas_TaxTable.to_csv(fnT, sep='\t')
    
    '''
    #Create Neufit string to run progam in terminal
    #Can now run Neufit entirely from this notebook so don't need this command
    py = 'python2 /Users/cguccione/Documents/Classes/Kit_Lab/Tools/neufit/neufit.py '
    fnD = finalFilename + '_data.csv'
    fnT = finalFilename + '_taxonomy.csv'
    neufit_command = py + fnD + ' -t ' + fnT
    
    return(neufit_command)
    Outputs: ![name]_data.csv = Creates a data.csv file to be input into Neufit located at 'neufit_input'
                ![name]_taxonomy.csv = Creates a taxonomy.csv file to be input into Neufit located at 'neufit_input'
    '''
    
    return(fnD, fnT)


def biom_addmetatax_custom_tcga_ehn(datasetName, biomFilename, finalFilename):
    '''Inputs: datasetName = filename of dataset ex. hutchKrakenAlex_biom -- found in set paths above
               biomFilename = filename of biom file found in datasetName path above
               finalFilename = filename to be used in _data.csv file for neuFit -- in this case will be normal or not
       Outputs: fnD_hn = File location of head_neck_[name]_data.csv
                fnD_e = File location of esophgous_[name]_data.csv
                !head_neck_[name]_data.csv = Creates a data.csv file to be input into Neufit located at 'neufit_input'
                !esophgous_[name]_data.csv = Creates a data.csv file to be input into Neufit located at 'neufit_input'
                fnT = File location of [name]_taxonomy.csv
                ![name]_taxonomy.csv = Creates a taxonomy.csv file to be input into Neufit located at 'neufit_input'
       * Imports biom file -> pandas dataframe -> splits apart hn from esoph -> 2 data.csv, 1 taxonomy.csv as desired by neufit
       * This custom because we need to seperate esophgous from head and neck patients from specific dataset'''

    #Make filename and import all data into pandas df
    fullFilename = datasetName + '/' + biomFilename
    featureTable = load_table(fullFilename) #Loads biom file: https://biom-format.org/documentation/generated/biom.load_table.html
    pandas_featureTable = pd.DataFrame(featureTable.matrix_data.toarray(), featureTable.ids('observation'), featureTable.ids())  #Turn biom file into pandas df
    meta = pd.read_csv(tcgaEhnWGSgreg_meta, sep = '\t') #Import the metadata into a pandas df
    pandas_TaxTable = pd.read_csv(tcgaEhnWGSgreg_taxa) #Import the taxonomy into a pandas df

    #Custom : Seperate the Head and Neck Samples from Esophgous Samples in the biom file
    ##In case Kit asks for it: biomDF has both hn and esoph samples combined
    hn_df = pd.DataFrame() #Biom file with just the head and neck samples
    e_df = pd.DataFrame() #Biom file with just the esophagus samples

    for sample in pandas_featureTable.columns.values: #Loop through all the samples in the biom file

        #Use the meta data to determine what cancer type the sample correlates to
        cancerType = meta[meta['sample_name'] == sample].reset_index().loc[0, 'primary_site'] #Finds the sample name which matches our i, then finds the corresponding primary site in that row

        #Sort each cancer type into approriate dataframes
        if cancerType == 'Head and Neck':
            hn_df[sample] = pandas_featureTable.loc[:, sample]
        elif cancerType == 'Esophagus':
            e_df[sample] = pandas_featureTable.loc[:, sample]
        else:#This should never print - just a saftey check
            print("something is off in the dataset -- not just Head and Neck and Esophagus")

    #Create _data.csv and _taxonomy.csv files for head & neck and esophgous

    #Head and Neck
    fnD_hn = neufit_input_path + '/' + 'head_neck_' + finalFilename + '_data.csv'
    fnT_hn = neufit_input_path + '/' +  'head_neck_' +  finalFilename + '_taxonomy.csv'
    pandas_TaxTable.to_csv(fnT_hn, sep='\t')
    pandas_featureTable.to_csv(fnD_hn, sep='\t')

    #Esophgous
    fnD_e = neufit_input_path + '/' + 'esophagus_' + finalFilename + '_data.csv'
    fnT_e = neufit_input_path + '/' +  'esophagus_' +  finalFilename + '_taxonomy.csv'
    pandas_TaxTable.to_csv(fnT_e, sep='\t')
    pandas_featureTable.to_csv(fnD_e, sep='\t')


    return(fnD_hn, fnD_e, fnT_hn, fnT_e)