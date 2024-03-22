# calculates rmsd using pymol
import sys
from pymol import cmd

# Load your proteins
cmd.load(sys.argv[1],"protein1")
cmd.fetch(sys.argv[2],"protein2")
cmd.align("protein1","protein2")

# Create selections for alpha carbons of functional sites
cmd.select("protein1_ca", "protein1 and i. 811-824+845+848+850-863+22+27-33+35+36+50-64+118-126+195+196+199+241-246+462+545-547+556-572+574-587+590+591+594+595+598+599+602+358+633+452 and name CA") 
cmd.select("protein2_ca", "((protein2 and polymer.protein within 10 of polymer.nucleic) and name CA) or (protein2 and i. 537+812+631 and n. CA)")

# Ensure the selections have the same number of atoms
if(cmd.count_atoms("protein1_ca") != cmd.count_atoms("protein2_ca")):
    print("Error: selections don't have the same number of atoms!")
else:
    # Get the model for each selection
    model1 = cmd.get_model("protein1_ca")
    model2 = cmd.get_model("protein2_ca")

    # Initialize a variable for the total RMSD
    total_rmsd = 0.0

    # Calculate per-residue RMSD
    for atom1 in model1.atom:
        for atom2 in model2.atom:
        # Calculate the squared distance between the atoms
            if int(atom1.resi)<699:
                if int(atom2.resi)==int(atom1.resi)+179:
                    dx = atom1.coord[0] - atom2.coord[0]
                    dy = atom1.coord[1] - atom2.coord[1]
                    dz = atom1.coord[2] - atom2.coord[2]
                    squared_distance = dx*dx + dy*dy + dz*dz
                    # Calculate the RMSD for this residue
                    rmsd = (squared_distance)**0.5
                    print(f"Residue {atom1.resi} and {atom2.resi}: RMSD = {rmsd}")
            if int(atom1.resi)>720:
                if int(atom2.resi)==int(atom1.resi)-720:
                    dx = atom1.coord[0] - atom2.coord[0]
                    dy = atom1.coord[1] - atom2.coord[1]
                    dz = atom1.coord[2] - atom2.coord[2]
                    squared_distance = dx*dx + dy*dy + dz*dz
                    # Calculate the RMSD for this residue
                    rmsd = (squared_distance)**0.5
                    print(f"Residue {atom1.resi} and {atom2.resi}: RMSD = {rmsd}")

        # Add this RMSD to the total
        total_rmsd += rmsd

    # Calculate the average RMSD
    average_rmsd = total_rmsd / cmd.count_atoms("protein1_ca")
    print(f"Average RMSD: {average_rmsd}")
