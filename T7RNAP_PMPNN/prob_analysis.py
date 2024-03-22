# This code plots the amino acid probabilities for generated sequences.

#import the necessary packages to read and plot data
import numpy as np
import plotly.express as px

address="/Volumes/lab-debenedictise/home/users/acarb/T7RNAP_PMPNN/230720_multi_pmpnn_4RNP_scores_probs_corrected/outputs/probs/"
fname="4RNP.npz"
scorefile = f'{address+fname}'#adjust as needed

#loading the scorefile, containng multiple numpy arrays
data = np.load(scorefile)

# Compute the size of the third part
part_size = len(data['probs'][1]) // len(data['chain_order'][0])

#plot amino acid scores as probabilities
alphabet = 'ACDEFGHIKLMNPQRSTVWYX' #amino acids 

fig = px.imshow(np.exp(data['probs'][:,:part_size,:]).mean(0).T, #for plotting the logarithmic probability change 'probs' to 'log_probs'
                labels=dict(x="positions", y="amino acids", color="probability"),
                y=list(alphabet),
                template="simple_white"
               ) # plot the probs matrix (conntaining the per residue score for the whole complex) and label the axes
fig.update_xaxes(side="top",
                dtick = 20)
fig.update_yaxes(tickmode='linear')

#save figure
output = f'{address+"probs.png"}' #adjust as needed
fig.write_image(output)
fig.show()
