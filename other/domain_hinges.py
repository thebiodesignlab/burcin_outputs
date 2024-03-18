from pymol import cmd, stored, math
import sys
sys.path.append('/Users/acarb/Dropbox (The Francis Crick)/DeBenedictisE/burcin/crick_biodesign/postdoc/gitacarbn')
import loadBfacts as lB
PDBid='1OMP'
sourcef="cross1.txt"
sourceh="hinges1.txt"
inFile = open(sourceh, 'r')
hinges=inFile.read()

hingesc = '+'.join(hinges.split())
print(hingesc)

cmd.fetch(PDBid)
cmd.remove('solvent')
lB.loadBfacts(PDBid,startaa=1,source=sourcef,visual="Y")
cmd.select("sele","resi "+hingesc+" and n. CA")
cmd.show("sphere","sele")
cmd.set("sphere_color","chartreuse")