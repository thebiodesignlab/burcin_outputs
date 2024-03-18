from pymol import cmd
import os
#NMRfiles=['1A03','2LGV','7UO6','7R0R','2KB1','7VIL']
NMRfiles=['8F4V','7ZE0','4TRX','1CIS','2ABD']

#os.mkdir("firstchains")
for n in NMRfiles:
    cmd.fetch(n)
    chs=cmd.get_chains(n)
    crm="+".join((chs[1:len(chs)]))
    if len(crm)>0:
        cmd.remove("chain "+crm)
    cmd.split_states(n)
    ss=cmd.get_object_list(n+'_*')
    for s in ss:
        cmd.save('./firstchains/'+s+'.pdb',s)
    del s
cmd.quit()
