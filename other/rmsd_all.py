# calculates rmsd using pymol
from pymol import cmd
import os

rmsds = []
PDB1='/Volumes/lab-debenedictise/home/users/acarb/230424_onlycpGFP/cpGFP/ranked_0.pdb'
cmd.load(PDB1)
cmd.set_name("ranked_0","refGFP")
all_files = os.listdir('/Volumes/lab-debenedictise/home/users/acarb/230317_MBP_EGFP_GPU_nomin/oldAF2/')
all_dir = [d for d in all_files if os.path.isdir('/Volumes/lab-debenedictise/home/users/acarb/230317_MBP_EGFP_GPU_nomin/oldAF2/%s'%d)]
num_dir = len(all_dir)
print(num_dir)

for i in range(num_dir):
    PDB2= '/Volumes/lab-debenedictise/home/users/acarb/230317_MBP_EGFP_GPU_nomin/oldAF2/fusion%d/fusion%d/ranked_0.pdb' % (i, i)
    cmd.load(PDB2)
    #cmd.align("refMBP", "ranked_0")
    rmsd=cmd.align("refGFP", "ranked_0", cutoff=2.0, cycles=5, gap=-10.0, extend=-0.5, max_gap=50, object=None, matrix='BLOSUM62', mobile_state=0, target_state=0, quiet=1, max_skip=0, transform=1, reset=0)
    rmsds.append(rmsd[3])
    cmd.delete("ranked_0")
    print(i,rmsd[0])
