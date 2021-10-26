# Caitlin Guccione, 08-25-2021

from matplotlib import pyplot
from neuevo.neuplot import neufit_plot

def custom_color_plot(occurr_freqs, n_reads, n_samples, r_square, beta_fit, file_header):
    '''Inputs: file_header = file path for neufit outpus: The path+dataGroup name + date stamp to be used for all parts of Neufit run: ex. /home/cguccion/NeutralEvolutionModeling/neufit_output/hutchKrakenAlex_combined_2021-08-26_13:24:22
               occurr_freqs = pandas df that Neufit created in the orginal program but didn't specifcally output orginally
                    - Headers of csv: otu_id, mean_abundance, occurrence, Kingdom, Phylum, Class, Order, Family, Genus, Species, predicted_occurrence, lower_conf_int, upper_conf_int
                    - I used this occur_freqs df in order to figure out which species is the most non-neutral
                    - I believe this is what Neufit uses to physical plot the dots on the graph as well
               n_reads = Number of reads (from orginal program)
               n_samples = Number of smaples (from orginal program)
               r_square = R^2 value (from orginal program)
               beta_fit = stats on the preformance of the model (from orginal program)
        Outputs: ! The plot that Neufit outputs + the colored dots - will print to the screen but not save when run
                 fn = default filename for neufit plot : _NeutralFitPlot_withNonNeutralColors.png
                 *Creates the plot for Neufit + with colored dots '''

    #Create orginal black plot
    neufit_plot(occurr_freqs, n_reads, n_samples, r_square, beta_fit, file_header)

    #If you want to add a new custom color, follow the steps below on what to add where (more info above as well)

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

    pyplot.plot(MA_P_Proteobacteria, O_P_Proteobacteria, 'o', markersize=6, fillstyle='full', color='green', label="Proteobacteria")
    pyplot.plot(MA_G_Streptococcus, O_G_Streptococcus, 'o', markersize=6, fillstyle='full', color='orange', label="Streptococcus")
    pyplot.plot(MA_P_Bacteroidetes, O_P_Bacteroidetes, 'o', markersize=6, fillstyle='full', color='purple', label="Bacteroidetes")
    pyplot.plot(MA_S_pylori, O_S_pylori, 'o', markersize=6, fillstyle='full', color='blue', label='H.pylori')


    plot_fn = file_header + '_NeutralFitPlot_withNonNeutralColors.png'

    return(plot_fn)

def save_plot(fn):
    ##Just saves the plot given a filename
    pyplot.legend(loc="center left")
    pyplot.savefig(fn)