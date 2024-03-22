# Instead of a 2D heatmap this script calculates for each generated sequence the average
# of sequence similarity with the rest of the sequences.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import squareform


from scipy.cluster.hierarchy import linkage, dendrogram, fcluster


# Step 1: Read and parse the data from the text file
address="/Volumes/lab-debenedictise/home/users/acarb/T7RNAP_PMPNN/singlePDB/230913_pmpnn_1CEZ/outputs/seqs/"

filename = address+'extract_data.txt'

xs, ys, values = [], [], []

with open(filename, 'r') as f:
    for line in f:
        if len(line)>4:
            print(line)
            print(len(line))
            x, y, value = line.split()
            xs.append(int(x))
            ys.append(int(y))
            xs.append(int(y))
            ys.append(int(x))
            values.append(float(value.rstrip('%')))
            values.append(float(value.rstrip('%')))

# Convert lists to numpy arrays
xs = np.array(xs)
ys = np.array(ys)
values = np.array(values)

# Create a grid for heatmap
grid = np.zeros((xs.max() + 1, ys.max() + 1))
for x, y, value in zip(xs, ys, values):
    grid[x, y] = value

# Step 2: Create the heatmap
#plt.imshow(grid, cmap='bwr', interpolation='nearest', origin='lower', vmin=0, vmax=50)
#plt.colorbar(label='Percentage')
#plt.xlabel('X-axis')
#plt.ylabel('Y-axis')
#plt.title('Heatmap from data')
#plt.show()

# Create a mask for the diagonal
mask = np.eye(grid.shape[0], grid.shape[1], dtype=bool)

# Create a masked array 
masked_grid = np.ma.array(grid, mask=mask)

# Compute the mean along the y-axis (axis=1) while ignoring masked values
averages = masked_grid.mean(axis=1).data

# Plot the averages
plt.plot(averages, marker='o', linestyle='-')
plt.xlabel('Sequence No')
plt.ylabel('Pairwise Similarity')
plt.title('Averages of Pairwise Similarity')
plt.grid(True)
plt.xticks(ticks=np.arange(0, grid.shape[0], 50))
plt.tight_layout()
plt.show()

# Set diagonal values to NaN
#for i in range(grid.shape[0]):
#    grid[i, i] = np.nan
# Plot each row on top of each other
#for i in range(grid.shape[0]):
#    plt.plot(grid[i, :], label=f'X={i}')

#plt.xlabel('Sequence No')
#plt.ylabel('Pairwise Similarity')
#plt.title('4RNP-Generated Sequences')
#plt.legend()
#plt.grid(True)
#plt.xticks(ticks=np.arange(0,grid.shape[1],50))  # Set x-ticks based on the number of Y values
#plt.tight_layout()
#plt.show()

# Create a clustered heatmap using seaborn
#sns.clustermap(grid, method='ward', cmap='coolwarm', figsize=(8, 6))

#plt.show()




