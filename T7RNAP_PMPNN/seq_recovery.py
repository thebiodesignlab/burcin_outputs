# This code calculates sequence recoveries for PMPNN generated sequences.
import re
import matplotlib.pyplot as plt
import numpy as np


address="/Volumes/lab-debenedictise/home/users/acarb/T7RNAP_PMPNN/singlePDB/230913_pmpnn_1CEZ/outputs/seqs/"
fname="1CEZ_A.fa"

# Open the file
with open(f'{address+fname}', 'r') as file:
    data = file.read()

# Regular expression to match the seq_recovery values
pattern = re.compile(r'seq_recovery=(\d+\.\d+)')
pattern2 = re.compile(r'global_score=(\d+\.\d+)')

# Find all matches in the file
matches = pattern.findall(data)
matches2 = pattern2.findall(data)

# Convert matches to floats
seq_recovery_values = [float(match) for match in matches]
global_score_values = [float(match2) for match2 in matches2]

# Print the matches
for i, value in enumerate(global_score_values):
    print(f'global_score_values for instance {i+1}: {value}')
"""
# Create an array for the x values (instance numbers)
x_values = range(1, len(seq_recovery_values) + 1)

# Create the figure and the first axes
fig, ax1 = plt.subplots(figsize=(10, 6))


# Plot the first array on the first axes
ax1.scatter(x_values,seq_recovery_values, color='blue',label='seq_recovery')
ax1.set_xlabel('Sequence Number',fontweight='bold')
ax1.set_ylabel('Sequence Recovery', color='blue',fontweight='bold')
ax1.tick_params(axis='y', labelcolor='blue')

# Create the second axes that shares the same x-axis as ax1
ax2 = ax1.twinx()

# Plot the second array on the second axes
ax2.scatter(x_values,global_score_values[1:], color='red',label='Global score')
ax2.set_ylabel('Global Score', color='red',fontweight='bold')
ax2.tick_params(axis='y', labelcolor='red')

# Add a title
fig.suptitle('4RNP Templated Sequences', fontweight='bold')

# Get the legend handles and labels from each axes
handles1, labels1 = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# Create a combined legend
plt.legend(handles1 + handles2, labels1 + labels2, loc="upper right")
plt.grid(True)
plt.show()"""

plt.figure(figsize=(6, 6))
plt.scatter(seq_recovery_values,global_score_values[1:])
#diag_values = np.linspace(0, max(max(x_values), max(seq_recovery_values)), 100)
#plt.plot(diag_values, diag_values, 'r--')
plt.grid(True)
plt.xlabel('Sequence recovery')
plt.ylabel('Global score')
plt.title(fname[0:4]+'-CA')
plt.show()

# Plot the seq_recovery values
#plt.figure(figsize=(6, 6))
#plt.scatter(x_values, seq_recovery_values)
#plt.title('Sequence recovery values')
#plt.xlabel('Instance')
#plt.ylabel('Sequence recovery')
#plt.grid(True)
#plt.show()
