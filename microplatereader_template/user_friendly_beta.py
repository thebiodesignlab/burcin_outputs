from helper_functions import main_make_graphs_wrapper, load_sample_names_to_legend
import pdb
import csv
import pickle
import random
            
number_of_plates = 1        # number of 96-well plates you used
dry_out_time = 50           # only include the first X timepoints
x_limits = (0, 0.6)         # x limits (OD) for the OD-lum graph
OD_lum_metric_cutoff = 0.3  # the OD at which you look at the luminescence across samples

    
# define mapping between plasmids, legend annotation, and color
# 'plasmid': ('legend annotation', color)
# see https://matplotlib.org/examples/color/named_colors.html
legend_dict = {
'pAB165c': ('PAB165C PED14xS-WT', 'black'),
}

# Auto-load and color samples
legend_dict = load_sample_names_to_legend(legend_dict)

# These will be automatically plotted on every graph
list_of_controls= ['PAB165C PED14xS-WT']

# For each graph, make a tuple that's (samples to plot, graph_title), and add it to the list of graphs
graphs_to_make = []
graphs_to_make.append(  (legend_dict.keys(), "All Samples"))

samples = legend_dict.keys()
for s in samples:
    graphs_to_make.append(  ([x for x in legend_dict.keys() if s in x], s))

# call the function to make graphs
main_make_graphs_wrapper(number_of_plates, dry_out_time, OD_lum_metric_cutoff, x_limits, legend_dict, list_of_controls, graphs_to_make)
