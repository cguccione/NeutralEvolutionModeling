import os
import click
from .__init__ import cli
from nevo.utils import (biom_addmetatax_custom_tcga_ehn,
                          non_neutral_outliers)
from nevo.neutral_fit import neufit


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
    required=False,
    default=True,
    help='TODO')
@click.option(
    '--colored-graph',
    required=False,
    default=True,
    help='TODO')
@click.option(
    '--non-neutral',
    required=False,
    default=True,
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
                      norm_graph : bool = True,
                      colored_graph : bool = True,
                      non_neutral : bool = True,
                      non_save : bool = False,
                      full_non_neutral : bool = False):
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
    
    '''
    # run within wrapper
    neufit(fnData, fnTaxonomy, output_filename,
           dataset_type, custom_filename,
           norm_graph = norm_graph, colored_graph = colored_graph, non_neutral = non_neutral,
           non_save = non_save,
           full_non_neutral = full_non_neutral)