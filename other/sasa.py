# solvent accessible surface area
import __main__
__main__.pymol_argv = ['pymol','-qc']
import pymol
from pymol import cmd, stored

cmd.set('dot_solvent', 1)
cmd.set('dot_density', 3)

cmd.fetch('1OMP')  # use the name of your pdb file
stored.residues = []
cmd.iterate('name ca', 'stored.residues.append(resi)')

sasa_per_residue = []
count=0
for i in stored.residues:
   # sasa_per_residue.append(cmd.get_area('resi %s' % i))
    sasa_per_residue.append(cmd.get_sasa_relative('resi %s' % i))
    print(i,sasa_per_residue[count])
    count+=1

print(sum(sasa_per_residue))
print(cmd.get_area('all'))  # just to check that the sum of sasa per residue equals the total area
