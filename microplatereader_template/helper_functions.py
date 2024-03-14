import pdb
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from scipy.optimize import curve_fit
from scipy.integrate import quad
import statistics
import matplotlib.patches as mpatches
import csv
from dateutil import parser
import sys
import os
import pickle
import random
import glob

# auto-load sample names into the legend dictionary, and assign them colors
def load_sample_names_to_legend(legend_dict):
    unique_sample_names = []
    with open('Manifest.csv', newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            row = [x for x in row if x != '']
            if len(row) > 2:
                if row[2] not in unique_sample_names:
                    unique_sample_names.append(row[2])

    # auto color samples
    all_colors = good_colors()
    for i in range(len(unique_sample_names)):
        x = unique_sample_names[i]
        if x == 'pAB165c':
            continue
        legend_dict[x] = (x, all_colors[i])
        
    # save legend dict to file to be read in by the analysis script
    pickle.dump(legend_dict, open( "legend_dict.p", "wb" ) )
    
    return legend_dict
        
# import colors
from matplotlib import colors as mcolors
def good_colors():
    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name) for name, color in colors.items())
    all_colors = []
    for (h,s,v), name in by_hsv:
        if s > .3 and v > .3:   # eliminate things too close to white
            all_colors.append(name)
    random.shuffle(all_colors)
    return all_colors


class Logger(object):
    # causes everything printed to stdout to also go to output.txt
    # https://stackoverflow.com/questions/14906764/how-to-redirect-stdout-to-both-file-and-console-with-scripting
    def __init__(self):
        self.terminal = sys.stdout
        if not os.path.exists("output.txt"):
            fn = open("output.txt", 'w+')
            fn.close()
        os.remove("output.txt")
        self.log = open("output.txt", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass    

sys.stdout = Logger()

INDIVIDUAL_WELL_GRAPHS = 0      ## generates time vs. OD and lum graphs
REGENERATE_GRAPHS = 1           ## generates OD vs. lum graphs

# disable open windows warning
plt.rcParams.update({'figure.max_open_warning': 0})

min_t = 0 
       
# calculate the lum at OD 0.3
# metric of suppression efficiency
def lum_at_od(OD, lum, target_OD):
    if max(OD) < target_OD:
        return None
    return np.interp(target_OD, OD, lum)

# OD at 6 hours
def growth_metric_2 (OD,t):
    time_delta = timedelta(hours=6)
    
    datetime_at_hr = t[0] + time_delta 
    float_datetime_at_hr = datetime_to_float(datetime_at_hr, min_t)
    
    # convert to float to interpolate
    float_t = [datetime_to_float(x, min_t) for x in t]
    
    OD_at_hr = np.interp(float_datetime_at_hr, float_t, OD)

    return OD_at_hr
    
# calculate time between OD 0.2 and OD 0.4
def growth_rate_metric(OD, t):
    if max(OD) < .4:
        return None

    # convert to float to interpolate
    float_t = [datetime_to_float(x, min_t) for x in t]
    
    time_2 = np.interp(0.2, OD, float_t)
    time_4 = np.interp(0.4, OD, float_t)
    
    # convert back to datetimes
    time_2_datetime = float_to_datetime(time_2, min_t)
    time_4_datetime = float_to_datetime(time_4, min_t)
    
    td = time_4_datetime - time_2_datetime
    
    return td.days * 24 * 60 + td.seconds / 60

def load_data(file):
    # loads data out of files that have timestamps across the top and well names down the columns
    data = {}
    with open(file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            well_num = ''
            for reading in row.items():
                (t, x) = reading
                if t == "Time":
                    data[x] = []
                    well_num = x
                    continue
                time = parser.parse(t)
                value = float(x)
                data[well_num].append((time, value))
    return data
    
def make_x_axis_readable(plt, time_range):
    # make there be fewer labels so that you can read them
    deltas = [t - time_range[0] for t in time_range]
    labels = [d.seconds/60/60 + d.days*24 for d in deltas]
    labels_sparse = [int(labels[x]) if x % 10 == 0 else '' for x in range(len(labels))]
    plt.xticks(time_range, labels_sparse)
    locs, labels = plt.xticks() 
    plt.xlabel('Hours')

def datetime_to_float(date, min_t):
    # minutes after the start of the experiment
    return (date - min_t).total_seconds()/60.
    
def float_to_datetime(fl, min_t):
    return min_t + timedelta(minutes = fl)

def perdelta(start, end, delta):
    curr = start
    while curr < end:
        yield curr
        curr += delta        

def plot_color(well, legend_dict):    
    for k in legend_dict.keys():
        if well in k:
            return legend_dict[k][1]  
    
    print("Error: trying to plot ", k, " but color cannot be found")

def auto_detect_file_name(n):
    '''Automatically determines the format of the .csv filenames'''
    extension = 'csv'
    result = glob.glob('*.{}'.format(extension))
    contains_string = 'reader_plate_' + str(n)
    my_files = [x for x in result if contains_string in x]
    abs = [x for x in my_files if 'abs' in x][0]
    lum = [x for x in my_files if 'lum' in x or 'supp_3_high' in x][0]
    return(abs, lum)

def load_file(plate_number, loaded_data):
    global min_t
    loaded_data[str(plate_number)] = {}
    
    abs_file, lum_file = auto_detect_file_name(plate_number)
    abs = load_data(abs_file)
    fluor = load_data(lum_file)

    #pdb.set_trace()
    
    discrete_downsample = {}
    calculated_values = {}
        
    for well in fluor.keys():
        
        # fetch values
        abs_t = [t for (t, y) in abs[well]]
        abs_y = [y for (t, y) in abs[well]]
        fluor_t = [t for (t, y) in fluor[well]]
        fluor_y = [y for (t, y) in fluor[well]]
        
        time_range = abs_t
        min_t = time_range[0]
        max_t = time_range[-1]
            
        # make plot
        if(INDIVIDUAL_WELL_GRAPHS):     
            # ensure that all subplots have the same x axis (time)
            i_times = [t for t in perdelta(min_t, max_t, timedelta(minutes=10))]
            i_times_float = [datetime_to_float(t, min_t) for t in i_times]
    
            # declare figure
            fig = plt.figure()
            #fig.suptitle(suppressor_name + ' - ' + well, fontsize=14)
        
            ## absorbance
            ax = plt.subplot(3, 1, 1)
            plt.xlim((min_t, max_t))
            plt.xticks([])
            plt.title( 'Plate 1 - ' + well, fontsize=14)
            plt.plot(abs_t, abs_y, 'r.')
            plt.ylabel('abs')
            
            ## fluorinescence
            ax = plt.subplot(3, 1, 2)
            plt.xlim((min_t, max_t))
            plt.xticks([])
            plt.plot(fluor_t, fluor_y, 'b.')
            plt.ylabel('fluor')

            ## fluor/abs
            ax = plt.subplot(3, 1, 3)
            plt.xlim((min_t, max_t))
            make_x_axis_readable(plt, time_range)
    
            i_fluor = np.interp(i_times_float, [datetime_to_float(t, min_t) for t in fluor_t], fluor_y)
            i_abs = np.interp(i_times_float, [datetime_to_float(t, min_t) for t in abs_t], abs_y)
    
            i_normfluor = i_fluor/i_abs
            plt.plot(i_times, i_fluor/i_abs, 'g.')
            plt.ylabel('fluor/abs')
        
            plt.tight_layout()
            plt.savefig('graphs/' + plate_name + '_' + well + '_Time_fluor.png', dpi = 300)

        od_values = np.interp([datetime_to_float(t, min_t) for t in fluor_t],
                              [datetime_to_float(t, min_t) for t in abs_t],
                               abs_y)
        loaded_data[str(plate_number)][well] = (od_values, fluor_y, fluor_t)
        
    return loaded_data

def main_make_graphs_wrapper(number_of_plates,
                        dry_out_time,
                        OD_lum_metric_cutoff,
                        x_limits,
                        legend_dict,
                        list_of_controls,
                        graphs_to_make):
    
    for i in ['OD vs lum', 't vs OD', 't vs lum']:
        main_make_graphs(number_of_plates,
                        dry_out_time,
                        OD_lum_metric_cutoff,
                        x_limits,
                        legend_dict,
                        list_of_controls,
                        graphs_to_make,
                        i)
                        
def main_make_graphs(number_of_plates,
                        dry_out_time,
                        OD_lum_metric_cutoff,
                        x_limits,
                        legend_dict,
                        list_of_controls,
                        graphs_to_make,
                        whats_on_the_axes):
    
    # load data from all files
    loaded_data = {}
    plate_numbers = [i for i in range(1, number_of_plates+1)]
    for p in plate_numbers:
        print("Loading plate " + str(p))
        loaded_data = load_file(p, loaded_data)

    metrics = {}
    
    # each iteration of this loop makes one graph
    for (plasmids_to_plot, graph_title) in graphs_to_make:

        fig = plt.figure()
        ax = plt.subplot(1, 1, 1)
    
        csvfile = open('Manifest.csv', 'r')
        reader = csv.reader(csvfile)
        for row in reader:
            (plate, well, contents, type) =  row[:4]
            plate_num = plate.split()[-1]

            if (contents in plasmids_to_plot or contents in list_of_controls):
                (od_values, fluor, t) = loaded_data[plate_num][well]
            
                # don't use all the data: the plate dried out
                od_values = od_values[:dry_out_time]
                lum = fluor[:dry_out_time]
                t = t[:dry_out_time]
                #print("Usable data for", plate, well, "is", t[0], "to", t[-1], datetime_to_float(t[-1],t[0])/60, "hr")

                if type == "Test suppressed":
                    ls = ':'
                else:
                    ls = '-'
            
            
                #pdb.set_trace()
                
                if whats_on_the_axes == "OD vs lum":
                    plt.plot(od_values, lum,
                             marker = 'o', linestyle = ls, linewidth=1, markersize=3, color = plot_color(contents, legend_dict))
                    plt.ylabel('luminescence (AU)')
                    plt.xlabel('OD')
                    plt.xlim(x_limits)
                elif whats_on_the_axes == "t vs OD":
                    plt.plot(t, od_values,
                             marker = 'o', linestyle = ls, linewidth=1, markersize=3, color = plot_color(contents, legend_dict))
                    plt.ylabel('OD')
                    plt.xlabel('t')
                    make_x_axis_readable(plt, t)
                elif whats_on_the_axes == "t vs lum":
                    plt.plot(t, lum,
                             marker = 'o', linestyle = ls, linewidth=1, markersize=3, color = plot_color(contents, legend_dict))
                    plt.ylabel('lum')
                    plt.xlabel('t')
                    make_x_axis_readable(plt, t)
                
                plt.title(graph_title, fontsize=14)
        
                # calculate the lum at OD 0.3
                lum_metric = lum_at_od(od_values, lum, OD_lum_metric_cutoff)
                ## accumulate the vlues here
                if (contents, type) in metrics:
                    if lum_metric not in metrics[(contents, type)]:
                        metrics[(contents, type)].append(lum_metric)
                else: 
                    metrics[(contents, type)] = [lum_metric]

                
        csvfile.close() 
        # write the OD to fluor to a file
        # collect the data for the legend

        print("Generating graph: ", graph_title)

        patchList = []
        for key in legend_dict.keys():
                if key in plasmids_to_plot or key in list_of_controls:
                    data_key = mpatches.Patch(color=legend_dict[key][1], label=legend_dict[key][0])
                    patchList.append(data_key)

        plt.legend(handles=patchList)#, bbox_to_anchor=(1.04,1), borderaxespad=0, loc='upper left')

        if whats_on_the_axes == "OD vs lum":
            ax.set_yscale('log')
            plt.autoscale(enable=True, axis='y')
        
        #plt.tight_layout()
        folder = 'graphs_' + whats_on_the_axes.replace(" ", '')
        if not os.path.exists(folder):
            os.makedirs(folder)
        plt.savefig(folder + "/" + graph_title + '.png', dpi = 300)

    ## Print summary information for plotting in Prism
    if whats_on_the_axes == "OD vs lum":
        print("\nLum at OD = ", OD_lum_metric_cutoff)    
        for k in metrics.keys():
            print(",".join(k), end = '')
            print(", ", end = '')
            my_list = metrics[k]
            # Print the formatted floats in a single line with comma-separated values (CSV) as separator
            # handling None values gracefully
            formatted_list = [("{:.2f}".format(num) if num is not None else "None") for num in my_list]
            print(','.join(formatted_list))

            
