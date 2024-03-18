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
            values.append(100-float(value.rstrip('%')))
            values.append(100-float(value.rstrip('%')))

# Add the y=x and color=100% data
unique_xs = set(xs)
for x in unique_xs:
    xs.append(x)
    ys.append(x)
    values.append(0.0)  # adding 100% value

# Convert lists to numpy arrays
xs = np.array(xs)
ys = np.array(ys)
values = np.array(values)

# Create a grid for heatmap
grid = np.zeros((xs.max() + 1, ys.max() + 1))
for x, y, value in zip(xs, ys, values):
    grid[x, y] = value
condensed_matrix = squareform(grid)
#print(condensed_matrix)  # Outputs: [1. 2. 3.]
Z= linkage(condensed_matrix, method='complete')  # or 'complete'

#plt.figure(figsize=(10, 7))#
#dendrogram(Z)
#plt.title("4RNP-Generated Sequences")
#plt.ylabel("Distances (100%-Similarity)")
#plt.show()

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111)
dendro_data = dendrogram(Z, ax=ax)
# Using the default threshold for coloring
# Using the default threshold for coloring
color_threshold = 0.7 * max(Z[:, 2])

clusters = []

for xs, ys in zip(dendro_data['icoord'], dendro_data['dcoord']):
    if ys[1] > color_threshold:  # Check if the height exceeds the threshold
        # Extract leaf indices from the x coordinates
        leaves = [int(x) for x in xs]
        clusters.append(leaves)



def create_lookup_table(matrix_size):
    lookup = {}
    index = 0
    for i in range(matrix_size):
        for j in range(i + 1, matrix_size):
            lookup[index] = (i, j)
            index += 1
    return lookup

# Assuming grid is your original distance matrix
lookup_table = create_lookup_table(grid.shape[0])

# Print the clusters
for i, cluster in enumerate(clusters, 1):
    print(f"{cluster}")


# Display the dendrogram
plt.figure(figsize=(10, 7))
dendrogram(Z)
plt.show()

# Step 2: Create the heatmap
#plt.imshow(grid, cmap='bwr', interpolation='nearest', origin='lower', vmin=0, vmax=50)
#plt.colorbar(label='Percentage')
#plt.xlabel('X-axis')
#plt.ylabel('Y-axis')
#plt.title('Heatmap from data')
#plt.show()

# Create a mask for the diagonal
#mask = np.eye(grid.shape[0], grid.shape[1], dtype=bool)

# Create a masked array 
#masked_grid = np.ma.array(grid, mask=mask)

# Compute the mean along the y-axis (axis=1) while ignoring masked values
#averages = masked_grid.mean(axis=1).data

# Plot the averages
#plt.plot(averages, marker='o', linestyle='-')
#plt.xlabel('Sequence No')
#plt.ylabel('Pairwise Similarity')
#plt.title('Averages of Pairwise Similarity')
#plt.grid(True)
#plt.xticks(ticks=np.arange(0, grid.shape[0], 50))
#plt.tight_layout()
#plt.show()

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




