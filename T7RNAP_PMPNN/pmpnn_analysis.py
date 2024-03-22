# This code plots the amino acid probabilities for generated sequences.
def clear_all():
    # Get a list of all variables in the global scope
    globals_to_clear = [var for var in globals().copy() if var[0] != "_"]
    # Get a list of all variables in the local scope
    locals_to_clear = [var for var in locals().copy() if var[0] != "_"]
    
    # Delete variables from the global scope
    for var in globals_to_clear:
        del globals()[var]
    
    # Delete variables from the local scope
    for var in locals_to_clear:
        del locals()[var]

clear_all()

import numpy as np
import matplotlib.pyplot as plt

# Load data from .npz file
address="/Volumes/lab-debenedictise/home/users/acarb/T7RNAP_PMPNN/230720_multi_pmpnn_4RNP_scores_probs_corrected/outputs/probs/"
fname="4RNP.npz"
data = np.load(address+fname)

#Print the keys in the file
print("Keys in the npz file: ", data.files)

# Access and print the arrays stored in the .npz file
for key in data.files:
    #print(f"{key}:\n{data[key]}")
    print(data[key].shape)
    #plt.plot(data[key], label=key)

# Get the list of keys
keys = list(data.files)
# Create the figure and the first axes
fig, ax1 = plt.subplots(figsize=(10, 6))
x_data = np.arange(len(data[keys[0]]))

# Plot the first array on the first axes
ax1.scatter(x_data,data[keys[0]], color='blue',label=keys[0])
ax1.set_xlabel('Sequence Number',fontweight='bold')
ax1.set_ylabel('Score', color='blue',fontweight='bold')
ax1.tick_params(axis='y', labelcolor='blue')

# Create the second axes that shares the same x-axis as ax1
ax2 = ax1.twinx()

# Plot the second array on the second axes
ax2.scatter(x_data,data[keys[1]], color='red',label=keys[1])
ax2.set_ylabel('Global Score', color='red',fontweight='bold')
ax2.tick_params(axis='y', labelcolor='red')


# Add a title
fig.suptitle('PMPNN Scores for 4RNP Templated Sequences', fontweight='bold')

# Get the legend handles and labels from each axes
handles1, labels1 = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# Create a combined legend
plt.legend(handles1 + handles2, labels1 + labels2, loc="upper right")



# Plot the data using matplotlib
#plt.figure(figsize=(6, 6))

# Plot each pair of consecutive arrays
#for i in range(0, len(keys) - 1, 1): # here we step by 2
#    x_data = data[keys[i]]
#    y_data = data[keys[i + 1]]
#    plt.scatter(x_data, y_data, label=f'{keys[i]} vs {keys[i + 1]}')
#    plt.xlabel(f'{keys[i]}',fontweight='bold')
#    plt.ylabel(f'{keys[i + 1]}',fontweight='bold')


# Add a title
#plt.title('PMPNN Scores for 4RNP Templated Sequences',fontweight='bold')
# Add a grid
plt.grid(True)
plt.show()

# Change print options
np.set_printoptions(threshold=np.inf)

# Now when you print your array, it won't be truncated
print(data['mask'])
