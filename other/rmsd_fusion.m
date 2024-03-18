clear all; clc; close all;
%% Preparing MBP component(PDB2)
PDB2=pdbread('/Volumes/lab-debenedictise/home/users/acarb/230424_onlyMBP/MBP/ranked_0.pdb');
PDB2.Sequence.NumOfResidues=PDB2.Model(1).Terminal.resSeq;
PDB2.Sequence.ChainID='A';
modelStruct2=PDB2.Model(1);
isChain2Atom = [modelStruct2.Atom.chainID]' == 'A'; % atoms in chain1
chain2Atoms = modelStruct2.Atom(isChain2Atom); % data for atoms in chain1
%=== retrieve sequence represented in the file (cannot always rely on pdbStruct1.Sequence(c1).Sequence)
chain2BkAtoms = chain2Atoms(strmatch('CA', {chain2Atoms.AtomName}, 'exact'));
r2=[chain2BkAtoms.resName];
k2=cellstr(reshape(r2,3,[])')';
l2=join(string(k2),' ');
m2=convertStringsToChars(l2);
PDB2.Sequence.ResidueNames=m2;

all_files = dir;
all_dir = all_files([all_files(:).isdir]);
num_dir = numel(all_dir) 

for i=0:num_dir
filename=sprintf('/Volumes/lab-debenedictise/home/users/acarb/230317_MBP_EGFP_GPU_nomin/oldAF2/fusion%d/fusion%d/ranked_0.pdb',i,i); 
%% Preparing PDB1
try
PDB1=pdbread(filename);
catch
disp(["cannot read ",filename])
end
PDB1.Sequence.NumOfResidues=PDB1.Model(1).Terminal.resSeq;
PDB1.Sequence.ChainID='A';
modelStruct1=PDB1.Model(1);
isChain1Atom = [modelStruct1.Atom.chainID]' == 'A'; % atoms in chain1
chain1Atoms = modelStruct1.Atom(isChain1Atom); % data for atoms in chain1
%=== retrieve sequence represented in the file (cannot always rely on pdbStruct1.Sequence(c1).Sequence)
chain1BkAtoms = chain1Atoms(strmatch('CA', {chain1Atoms.AtomName}, 'exact')); 
r=[chain1BkAtoms.resName];
k=cellstr(reshape(r,3,[])')';
l=join(string(k),' ');
m=convertStringsToChars(l);
PDB1.Sequence.ResidueNames=m;

%% distance & RMSD
[Dist(i+1), RMSD(i+1)] = pdbsuperpose(PDB2,PDB1,'SeqAlign',false); %'Display',false
clear PDB1 modelStruct1 isChain1Atom chain1Atoms r k l m
end
save('dist2_rmsd','Dist','RMSD') 
disp('done')
%quit
