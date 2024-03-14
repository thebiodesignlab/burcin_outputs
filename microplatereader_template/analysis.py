import pickle
import pdb
import numpy as np
from functools import reduce

# OD as a function of time
grm = pickle.load( open( "growth_rate_metrics.p", "rb" ) )
for plasmid, type in grm:
    my_metrics = [x for x in grm[(plasmid, type)] if x is not None]        
    print(plasmid, type, np.mean(my_metrics), np.std(my_metrics), len(my_metrics))

# luminescence as a function of OD
metrics = pickle.load( open( "metrics.p", "rb" ) )
vals = {}
all_vals = {}

for i in range(5,55):
    for plasmid, type in metrics[i]:
        my_metrics = [x for x in metrics[i][(plasmid, type)] if x is not None]        
        vals[plasmid, type, i] = (np.mean(my_metrics), np.std(my_metrics), len(my_metrics))
        all_vals[plasmid,type,i] = my_metrics

def error_propagate_subtraction(a, b):
    ''' Returns (mean, variance) tuple of a-b'''
    (a_mean, a_var) = a
    (b_mean, b_var) = b
    c_mean = a_mean - b_mean
    c_var = np.sqrt(np.power(a_var,2) + np.power(b_var,2))
    return (c_mean, c_var)

def error_propagate_division(n, d):
    ''' Returns (mean, variance) tuple of n/d
    Variables named for numerator/demoninator'''
    (n_mean, n_var) = n
    (d_mean, d_var) = d
    c_mean = float(n_mean)/float(d_mean)
    c_var = np.sqrt(1/(np.power(d_mean,4)) * np.power(d_var,2)
                  + 1/(np.power(d_mean,2)) * np.power(n_var,2) )
    return (c_mean, c_var)

def combine_std (N1, N2):
    (X1, s1, n1) = N1
    (X2, s2, n2) = N2
    new_mean = (n1*X1+n2*X2)/(n1+n2)
    d1 = X1-new_mean
    d2 = X2-new_mean
    new_std = np.sqrt((n1*(s1*s1+d1*d1)+n2*(s2*s2+d2*d2))/(n1+n2))
    return (new_mean, new_std, n1+n2)

# sanity check tests  
#print(error_propagate_subtraction((243078.7721, 3869.782666),(1779.702961, 278.0639858)))
#print(error_propagate_division((1538.452605, 388.5257109),(241299.0691, 3879.759974)))

def calculate_percent_WT(test_induced, test_suppressed, wt):
    s_test = error_propagate_subtraction(test_induced, test_suppressed)
    s_wt =   error_propagate_subtraction(wt, test_suppressed)
    
    return error_propagate_division(s_test, s_wt)
        
def percent_WT_wrapper(sample_name, control_name, od_cutoff):
    test = vals[(sample_name, 'Test induced', od_cutoff)]
    supp = vals[(sample_name, 'Test suppressed', od_cutoff)]
    control = vals[(control_name, 'Control', od_cutoff)]
    v = calculate_percent_WT((test[0], test[1]), (supp[0], supp[1]), (control[0], control[1]))
    ## ADDED!
    print("Percent WT", sample_name, ",", v[0]*100, ",", v[1]*100,",", end = '')

    all_v_induced = all_vals[sample_name, 'Test induced', od_cutoff]
    all_v_control = all_vals[sample_name, 'Test suppressed', od_cutoff]
    for i in range(len(all_v_induced)):
        try:
            print(all_v_induced[i] - all_v_control[i],",", end = '')
        except:
            pass
    print(" ")


def percent_WT_of_control_wrapper(sample_name, control_name, od_cutoff):    
    n = vals[(sample_name, 'Control', od_cutoff)]
    d = vals[(control_name, 'Control', od_cutoff)]
    v = error_propagate_division((n[0], n[1]), (d[0], d[1]))
    print("Percent WT", sample_name, ",", v[0]*100, ",", v[1]*100)

def percent_control_wrapper(sample_name, control, od_cutoff):
    test = vals[(sample_name, 'Test induced', od_cutoff)]
    supp = vals[(sample_name, 'Test suppressed', od_cutoff)]
    v = calculate_percent_WT((test[0], test[1]), (supp[0], supp[1]), (control[0], control[1]))
    print("Percent WT", sample_name, ",", v[0]*100, ",", v[1]*100)
    #print(v[0], v[1], end=' ')
    #print(od_cutoff, v[0], v[1])
    
# load in legend dict that was used by user_friendly_beta.py
legend_dict = pickle.load( open( "legend_dict.p", "rb" ) )
del legend_dict["pAB165c"]  # remove the positive control

    
#for od_cutoff in range(5, 55):
if(1):
    od_cutoff = 30
    
    print('pAB165c,', vals[('pAB165c', 'Control', od_cutoff)])
    for i in [x for x in legend_dict.keys()]:# if 'CAGG' in x]:
        if '1e0' in i:
            continue
        percent_WT_wrapper(i, 'pAB165c', od_cutoff)

