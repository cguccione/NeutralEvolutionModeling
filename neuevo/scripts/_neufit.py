import os
import click
from .__init__ import cli
from neuevo.utils import (biom_addmetatax_custom_tcga_ehn,
                          non_neutral_outliers)
from neuevo.neufit import neufit


@cli.command(name='neufit')
@click.option(
    '--fnData',
    required=True,
    help='TODO')
@click.option(
    '--fnTaxonomy',
    required=True,
    help='TODO')
@click.option(
    '--output-filename',
    required=True,
    help='TODO')
@click.option(
    '--dataset-type',
    required=True,
    help='TODO')
@click.option(
    '--custom-filename',
    required=True,
    help='TODO')
@click.option(
    '--norm-graph',
    required=True,
    help='TODO')
@click.option(
    '--colored-graph',
    required=True,
    help='TODO')
@click.option(
    '--non-neutral',
    required=True,
    help='TODO')
@click.option(
    '--non-save',
    required=False,
    default=False,
    help='TODO')
@click.option(
    '--full-non-neutral',
    required=False,
    default=False,
    help='TODO')
def standalone_neufit(fnData : str,
                      fnTaxonomy : str,
                      output_filename : str,
                      dataset_type: str,
                      custom_filename : str,
                      norm_graph : bool,
                      colored_graph : bool,
                      non_neutral : bool,
                      non_save : bool = False,
                      full_non_neutral : bool = False):
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
    # run within wrapper
    neufit(fnData, fnTaxonomy, output_filename,
           dataset_type, custom_filename,
           norm_graph, colored_graph, non_neutral,
           non_save = non_save,
           full_non_neutral = full_non_neutral)