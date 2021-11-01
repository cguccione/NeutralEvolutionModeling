from matplotlib import pyplot
from nevo.neutral_fit_plot import neufit_plot

def custom_color_plot(occurr_freqs, n_reads, n_samples, r_square, beta_fit, file_header):
    '''Adds species/phylum specific coloring to the neutral evolution plot
    
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
    n_reads: int
        Total number of reads.
    n_samples: int
        Total number of samples.
    r_square: float
        R^2 value of the fit of data to neutral curve.
    beta_fit: lmfit.model.ModelResult object
        Holds the stats on the preformance of the model.
    
    Returns
    -------
    png
        A plot showing the the input data and how well it fits the neutral
        model with specific speices or phylums labled in different colors
    fn: str, path
        A filepath where the plot is stored. 
    
    Notes
    -----
    If you want to add another custom color:
        1. Create a new section with bacteria and Mean Abundance 
        and Occurance lists
        2. Add to the for loop at bottom of section as well
        3. Add custum bacteria to plot
    
    TODO
    ----
    - Add a better definition for returns, is there a more clear 
    way to describe the plot
    - Confirm the filename outputed is acutally used in code
    '''
    
    #Create orginal black plot
    neufit_plot(occurr_freqs, n_reads, n_samples, r_square,
                beta_fit, file_header)
    
    #If you want to add a new custom color, follow the steps below on 
    #what to add where (more info above as well)
    
    # --------- STEP 1 ---------
    # Create mean abundance and occurance list for bacteria 

    #Phylori
    MA_S_pylori = [] #Mean_Abundance, Species , phylori
    O_S_pylori = [] #Occurrance

    #Proteobacteria
    MA_P_Proteobacteria = []
    O_P_Proteobacteria = []

    #Streptococcus
    MA_G_Streptococcus = []
    O_G_Streptococcus = []

    #Bacteroidetes
    MA_P_Bacteroidetes = []
    O_P_Bacteroidetes = []

    # --------- STEP 2 ---------
    # Add to customized lists when hitting bacteria in neutral list 

    #for i, j in neutral.iterrows():
    for i,j in occurr_freqs.iterrows():
        if j['Species'] == 's__pylori':
            MA_S_pylori.append(j['mean_abundance'])
            O_S_pylori.append(j['occurrence'])
        if j['Phylum'] == 'p__Proteobacteria':
            MA_P_Proteobacteria.append(j['mean_abundance'])
            O_P_Proteobacteria.append(j['occurrence'])
        if j['Genus'] == 'g__Streptococcus':
            MA_G_Streptococcus.append(j['mean_abundance'])
            O_G_Streptococcus.append(j['occurrence'])
        if j['Phylum'] == 'p__Bacteroidetes':
            MA_P_Bacteroidetes.append(j['mean_abundance'])
            O_P_Bacteroidetes.append(j['occurrence'])

    # --------- STEP 3 ---------
    # Override coloring of dots in current plot with new colors

    pyplot.plot(MA_P_Proteobacteria, O_P_Proteobacteria, 'o', 
                markersize=6, fillstyle='full', color='green', 
                label="Proteobacteria")
    pyplot.plot(MA_G_Streptococcus, O_G_Streptococcus, 'o', 
                markersize=6, fillstyle='full', color='orange', 
                label="Streptococcus")
    pyplot.plot(MA_P_Bacteroidetes, O_P_Bacteroidetes, 'o', 
                markersize=6, fillstyle='full', color='purple', 
                label="Bacteroidetes")
    pyplot.plot(MA_S_pylori, O_S_pylori, 'o', markersize=6, 
                fillstyle='full', color='blue', label='H.pylori')
        
    plot_fn = file_header + '_NeutralFitPlot_withNonNeutralColors.png'
    
    return(plot_fn)

def save_plot(fn):
    '''Saves the plot given filename input
    
    Written by: Caitlin Guccione, 08-25-2021
    '''
    pyplot.legend(loc="center left")
    pyplot.savefig(fn)