from scipy.spatial import distance_matrix
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import math

address="/Volumes/lab-debenedictise/home/users/acarb/T7RNAP_PMPNN/230720_multi_pmpnn_4RNP_scores_probs_corrected/outputs/seqs/"
fname="blastp.out"

data = []

# Open the BLASTP output file
with open(f'{address+fname}', 'r') as file:
    # Read all lines into a list
    lines = file.readlines()

    # Count the lines
    sl = len(lines)
    print(f'{sl}')
    
    # Initialize an empty list to hold the data
    data = []

    # Iterate over the lines
    for line in lines:
        # Split the line into fields
        fields = line.strip().split('\t')
        # Extract the necessary data
        pident = float(fields[2]) / 100
        
        # Add the data to the list
        data.append(1-pident)


pident_array = np.array(data)
"""distance_matrix = pident_array.reshape((1001, 1001))


# Now, distance_matrix can be used for clustering
sns.heatmap(distance_matrix, cmap='viridis')
# Show the plot
plt.show()"""