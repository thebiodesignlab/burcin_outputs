from pymol import cmd, stored, math
import sys
sys.path.append('/Users/acarb/Dropbox (The Francis Crick)/DeBenedictisE/burcin/crick_biodesign/postdoc/gitacarbn')
import spectrumany as sp

	
def loadBfacts (mol,startaa=1,source="newBfactors.txt", visual="Y"):
	"""
	Replaces B-factors with a list of values contained in a plain txt file
	
	usage: loadBfacts mol, [startaa, [source, [visual]]]
 
	mol = any object selection (within one single object though)
	startaa = number of first amino acid in 'new B-factors' file (default=1)
	source = name of the file containing new B-factor values (default=newBfactors.txt)
	visual = redraws structure as cartoon_putty and displays bar with min/max values (default=Y)
 
	example: loadBfacts 1LVM and chain A
	"""
	obj=cmd.get_object_list(mol)[0]
	cmd.alter(mol,"b=-1.0")
	inFile = open(source, 'r')
	counter=int(startaa)
	bfacts=[]
	for line in inFile.readlines():	
		bfact=float(line)
		bfacts.append(bfact)
		cmd.alter("%s and resi %s and n. CA"%(mol,counter), "b=%s"%bfact)
		counter=counter+1
	if visual=="Y":
		cmd.show_as("cartoon",mol)
		#cmd.cartoon("putty", mol)
		#cmd.set("cartoon_putty_scale_min", min(bfacts),obj)
		#cmd.set("cartoon_putty_scale_max", max(bfacts),obj)
		#cmd.set("cartoon_putty_transform", 0,obj)
		#cmd.set("cartoon_putty_radius", 0.2,obj)
		#cmd.spectrum("b","rainbow", "%s and n. CA " %mol,minimum=0, maximum=0.1)
		sp.spectrumany("b","purple white green", "(all)",minimum=-5, maximum=5) #color structure
		#cmd.ramp_new("count", obj, [-5, 5], 'rainbow') #color colorbar
		cmd.ramp_new("count", obj, [-5, 0, 5], ['purple','white','green']) #color the colorbar accordingly
		cmd.recolor()

cmd.extend("loadBfacts", loadBfacts);