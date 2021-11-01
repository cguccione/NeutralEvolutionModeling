import pandas as pd
from biom import load_table 
import os
import numpy as np

def biom2data_tax(datasetName, biomFilename, finalFilename):
    '''Imports biom file -> pandas dataframe -> data.csv, taxonomy.csv 
    
    Written by: Caitlin Guccione, 08-25-2021
    
    Parameters
    ----------
    datasetName: str
        The filename of dataset ex. hutchKrakenAlex_biom 
    biomFilename: str
        The filename of biom file
    finalFilename: str
        The filename to be used as the header in _data.csv and
        _taxonomy.csv files
        
     Returns
     -------
    fnD: str
        The file location of [name]_data.csv
    csv file
        Creates[name]_data.csv, a data.csv file needed to run Neufit 
    fnT: str
        The file location of [name]_taxonomy.csv
    csv file
        Creates [name]_taxonomy.csv, a taxonomy.csv file needed
        to run Neufit
      
     TODO
     ----
     - Insure there isn't a better way to inpur filename
     - Make sure the path direction is clear for all inputs
     - Add more clear description of csv files 
     
      '''
    neufit_input_path = '/home/cguccion/NeutralEvolutionModeling/ipynb/data_tax_csv' 
    #location of _data.csv and _tax.csv files for Neufit input

    
    #Make filename and import data
    fullFilename = datasetName + '/' + biomFilename
    featureTable = load_table(fullFilename) 
    #https://biom-format.org/documentation/generated/biom.load_table.html
    
    #Create _data.csv
    pandas_featureTable = pd.DataFrame(featureTable.matrix_data.toarray(),
                                       featureTable.ids('observation'), 
                                       featureTable.ids())
    fnD = neufit_input_path + '/' + finalFilename + '_data.csv'
    pandas_featureTable.to_csv(fnD, sep='\t')
    
    #Create _taxonomy.csv
    pandas_TaxTable = pd.DataFrame(featureTable.metadata_to_dataframe('observation'))
    pandas_TaxTable.set_axis(['Kingdom', 'Phylum', 'Class', 'Order',
                              'Family', 'Genus', 'Species'], axis=1)
    fnT = neufit_input_path + '/' + finalFilename + '_taxonomy.csv'
    pandas_TaxTable.to_csv(fnT, sep='\t')
    
    return(fnD, fnT)

def biom_addMetaTax_customTCGAehn(datasetName, biomFilename, finalFilename):
    ''' Imports biom file -> pandas dataframe -> 
        splits apart head and neck from esoph -> 2 data.csv, 1 taxonomy.csv 
       * This custom because we need to seperate esophgous 
       from head and neck patients from specific dataset *
    
    Written by: Caitlin Guccione, 08-25-2021
    
    Parameters
    ----------
    datasetName: str
        The filename of dataset ex. hutchKrakenAlex_biom 
    biomFilename: str
        The filename of biom file
    finalFilename: str
        The filename to be used as the header in _data.csv and
        _taxonomy.csv files
        
     Returns
     -------
    fnD_hn: str
        The file location of head_neck_[name]_data.csv
    fnD_e: str
        The file location of esophagous_[name]_data.csv
    csv file
        Creates head_neck_[name]_data.csv, a data.csv file 
    csv file
        Creates esophagous_[name]_data.csv, a data.csv file 
    fnT: str
        The file location of [name]_taxonomy.csv
    csv file
        Creates [name]_taxonomy.csv, a taxonomy.csv file needed
        to run Neufit
      
     TODO
     ----
     - Experinment if there is a more efficent wayt to run this
         so that you don't have to create an custom function as 
         done so here. Consider finding ways to split the data in 
         a seperate pre-processing script
     - Add more clear description of csv files 
     - Check that all ouputs match the return values here
     
      '''
    
    #Define filenames for meta and taxa files
    tcgaEhnWGSgreg_ = '/home/cguccion/rawData/April2021_Greg_TCGA_WGS/raw_from_Greg'#Location of raw data
    tcgaEhnWGSgreg_meta = str('/home/cguccion/rawData/April2021_Greg_TCGA_WGS/meta_expansion/13722_20210405-101126-TCGA-WGS-Qiita-sample-metadata_esoph_hnc_metaExpand.txt')
    tcgaEhnWGSgreg_taxa = str(tcgaEhnWGSgreg_ + '/' + 'wol_gotu_taxonomy.csv')
    
    
    #Make filename and import all data into pandas df 
    fullFilename = datasetName + '/' + biomFilename
    featureTable = load_table(fullFilename) 
    #Loads biom file: https://biom-format.org/documentation/generated/biom.load_table.html
    pandas_featureTable = pd.DataFrame(featureTable.matrix_data.toarray(),
                                       featureTable.ids('observation'),
                                       featureTable.ids())  
                                        #Turn biom file into pandas df 
    meta = pd.read_csv(tcgaEhnWGSgreg_meta, sep = '\t') 
    #Import the metadata into a pandas df
    pandas_TaxTable = pd.read_csv(tcgaEhnWGSgreg_taxa) 
    #Import the taxonomy into a pandas df
    
    #Custom : Seperate the Head and Neck Samples from Esophgous 
    #Samples in the biom file
    hn_df = pd.DataFrame() #Biom file with just the head and neck samples
    e_df = pd.DataFrame() #Biom file with just the esophagus samples
    all_scc_df = pd.DataFrame() #Biom file with Squamous Cell carcinoma esophagus samples
    eac_df = pd.DataFrame() #Biom file with EAC samples
    
    for sample in pandas_featureTable.columns.values: 
        
        #Use the meta data to determine what cancer type the sample correlates to 
        cancerType = meta[meta['sample_name'] == sample].reset_index().loc[0, 'primary_site'] 
        #Finds the sample name which matches our i, 
        #then finds the corresponding primary site in that row
        
        #Sort each cancer type into approriate dataframes
        if cancerType == 'Head and Neck':
            hn_df[sample] = pandas_featureTable.loc[:, sample]
        elif cancerType == 'Esophagus':
            #All esoph cancers
            e_df[sample] = pandas_featureTable.loc[:, sample]
            
            #Break esoph cancer into Squamous cell carcinoma, NOS, 
            #Squamous cell carcinoma, keratinizing, NOS & Adenocarcinoma, 
            #NOS -- there are a few other ones in there than you have ignored 
            esophCancerType = meta[meta['sample_name'] == sample].reset_index().loc[0, 'primary_diagnosis'] 
            #Finds the sample name which matches our i, then 
            #finds the corresponding primary site in that row
            
            if esophCancerType == 'Squamous cell carcinoma, NOS' or 'Squamous cell carcinoma, keratinizing, NOS':
                all_scc_df[sample] = pandas_featureTable.loc[:, sample]
            elif esophCancerType == 'Adenocarcinoma, NOS':
                eac_df[sample] = pandas_featureTable.loc[:, sample]
            
        else:#This should never print - just a saftey check 
            print("something is off in the dataset -- not just Head and Neck and Esophagus")
        
    #Create _data.csv and _taxonomy.csv files for head & neck and esophgous
    
    neufit_input_path = '/home/cguccion/NeutralEvolutionModeling/ipynb/data_tax_csv' 
    #location of _data.csv and _tax.csv files for Neufit input
    
    #Head and Neck
    fnD_hn = neufit_input_path + '/' + 'head_neck_' + \
        finalFilename + '_data.csv'
    fnT_hn = neufit_input_path + '/' +  'head_neck_' + \
        finalFilename + '_taxonomy.csv'
    pandas_TaxTable.to_csv(fnT_hn, sep='\t')
    hn_df.to_csv(fnD_hn, sep='\t')
    
    #Esophgous
    fnD_e = neufit_input_path + '/' + 'esophagus_' + \
        finalFilename + '_data.csv'
    fnT_e = neufit_input_path + '/' +  'esophagus_' + \
        finalFilename + '_taxonomy.csv'
    pandas_TaxTable.to_csv(fnT_e, sep='\t')
    e_df.to_csv(fnD_e, sep='\t')
    
    ##Squamous Cell Carcinoma
    fnD_e_scc = neufit_input_path + '/' + 'esophagus_SCC_' + \
        finalFilename + '_data.csv'
    fnT_e_scc = neufit_input_path + '/' +  'esophagus_SCC_' + \
        finalFilename + '_taxonomy.csv'
    pandas_TaxTable.to_csv(fnT_e_scc, sep='\t')
    all_scc_df.to_csv(fnD_e_scc, sep='\t')
    
    ##EAC
    fnD_e_eac = neufit_input_path + '/' + 'esophagus_EAC_' + \
        finalFilename + '_data.csv'
    fnT_e_eac = neufit_input_path + '/' +  'esophagus_EAC_' + \
        finalFilename + '_taxonomy.csv'
    pandas_TaxTable.to_csv(fnT_e_eac, sep='\t')
    eac_df.to_csv(fnD_e_eac, sep='\t')
    
    
    return(fnD_hn, fnD_e, fnT_hn, fnT_e, fnD_e_scc, fnT_e_scc, fnD_e_eac, fnT_e_eac)

def non_neutral_outliers(file_header, occurr_freqs, dataset_type, non_save, threshold = 0.5):
    ''' Creates the most Non-neutral csv file 
    
    Written by: Caitlin Guccione, 08-25-2021
    
    Parameters
    ----------
    file_header: str, path
        Filepath for all nevo outputs. Includes path, data nickname and 
        time stamp.
    occurr_freqs: pandas df
        Df header: otu_id, mean_abundance, occurrence, Kingdom, Phylum, 
        Class, Order, Family, Genus, Species, predicted_occurrence, 
        lower_conf_int, upper_conf_int
    threshold: int, optional
        Autoset to 0.5, but determines which bacteria are considered 
        non-neutral strictly for this csv file
        
     Returns
     -------
     csv file
         [name]_NonNeutralOutliers.csv, holds all the species 
         furthest off the neutral curve
    
    TODO
    ----
    - See if there is a way to eliminate the Domain vs Kingdom 
    differences more cleanly 
    
    '''
    
    #Create dataframe
    standoutMicrobes = pd.DataFrame(columns = ('Difference off Neutral Model',
                                               'Kingdom', 'Phylum', 'Class',
                                               'Order', 'Family', 'Genus',
                                               'Species'))
    
    #Loop and find most non-neutral microbes based upon threshold
    row_count = 0
    for i,j in occurr_freqs.iterrows():
        diff = abs(j['occurrence'] - j['predicted_occurrence'])
        if diff > threshold:
            if dataset_type == 'TCGA_WGS':
                standoutMicrobes.loc[row_count] = [diff, j['Domain'],
                                                   j['Phylum'], j['Class'],
                                                   j['Order'], j['Family'],
                                                   j['Genus'], j['Species']]
            else:
                standoutMicrobes.loc[row_count] = [diff, j['Kingdom'],
                                                   j['Phylum'], j['Class'],
                                                   j['Order'], j['Family'],
                                                   j['Genus'], j['Species']]
            row_count +=1
    standoutMicrobes = standoutMicrobes.sort_values(by =['Difference off Neutral Model'],
                                                    ascending=False)
    
    #Display and export non-neutral microbes as csv
    print("\nTop NonNeutral Microbes")
    display(standoutMicrobes)
    print('=========================================================\n')
    
    fn = file_header + '_NonNeutral_Outliers.csv'
    
    if non_save == False:
        standoutMicrobes.to_csv(fn, sep = '\t')