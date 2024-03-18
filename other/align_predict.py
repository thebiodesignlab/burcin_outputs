from pymol import cmd, stored, util
import sys, os, glob, math
sys.path.append('/Users/acarb/Dropbox (The Francis Crick)/DeBenedictisE/burcin/crick_biodesign/postdoc/gitacarbn')
import radiusOfGyration as rg
import center_of_mass as cm

cmd.set('sphere_transparency', 0.5)
cmd.fetch('4O4B')

r=rg.rgyrate('4O4B')
cm.com('polymer', object='com', vdw=r)
print(r)
util.cbc()
cmd.disable('com')


os.chdir('/Users/acarb/Dropbox (The Francis Crick)/DeBenedictisE/burcin/crick_biodesign/postdoc/')
filename='test_e62db.result/*.pdb'
fl=sorted(glob.glob(filename))
count=1
for i in fl:
    cmd.load(i)
    j=i.split('/')
    r=rg.rgyrate(j[1][:-4])
    
    cmd.alignto('4O4B and chain B')
    cm.com('polymer', object='com'+str(count), vdw=r)
    print(r)
    util.cbc()
    cmd.disable('com'+str(count))

    count=count+1
  #  cmd.align(i,'4O4B and chain A')