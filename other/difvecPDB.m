%% This code compares GNM modes with the dif vector between two states of a protein
% and plot the overlaps.
%%
function [resnum,ndvs,dvs,MSF1f,MSF2f]=difvecPDB(fname1,chain1,fname2,chain2)
% fnames 1 and 2 are PDB codes
% and chains 1 and 2 are the corresponding protein chains
f1=getpdb(fname1);
f2=getpdb(fname2);
[Dist, RMSD, Transf, PDB2TX] = pdbsuperpose(f1, f2);
clear f2
f2=PDB2TX;
atom1=length(f1.Model.Atom);
atom2=length(f2.Model.Atom);
count=0;
count2=0;

[MSF1i,MSF1a,resnum]=GNM(fname1,[1:10],chain1);
[MSF2i,MSF2a,resnum2]=GNM(fname2,[1:10],chain2);

for i=1:atom1
    if contains(chain1,f1.Model.Atom(i).chainID)==1 && ...
            strcmp('CA',f1.Model(1).Atom(i).AtomName)
            count=count+1;
            x(count,1)=f1.Model(1).Atom(i).X;
            y(count,1)=f1.Model(1).Atom(i).Y;
            z(count,1)=f1.Model(1).Atom(i).Z;
            sq=f1.Model(1).Atom(i).resName;
            sq1(count,1)=aminolookup(sq);
    end
end
for i=1:atom2
    if contains(chain2,f2.Model.Atom(i).chainID)==1 && ...
            strcmp('CA',f2.Model(1).Atom(i).AtomName)
            count2=count2+1;
            x(count2,2)=f2.Model(1).Atom(i).X;
            y(count2,2)=f2.Model(1).Atom(i).Y;
            z(count2,2)=f2.Model(1).Atom(i).Z;
            sq2=f2.Model(1).Atom(i).resName;
            sq12(count2,1)=aminolookup(sq2);
    end
end

if count==count2
    X=x(:,2)-x(:,1);
    Y=y(:,2)-y(:,1);
    Z=z(:,2)-z(:,1);
    dvs=sqrt(X.^2+Y.^2+Z.^2);
    MSF1f=MSF1i;
    MSF2f=MSF2i;
else
    disp('Please match the residue numbers of two proteins!')
    rdn=[strcat(sq1)];
    irdn=[strcat(sq12)];
    [Score, Alignment] =nwalign(rdn',irdn')
    xa=zeros(length(Alignment),2);
    for k=1:10
        MSF1f(k,:)=zeros(1,length(Alignment));
        MSF2f(k,:)=zeros(1,length(Alignment));
    end
    ya=xa; za=xa;
    for j=1:2
        count3=0;
    
        for i=1:length(Alignment)
           if strcmp(Alignment(2*j-1,i),'-')==0
               count3=count3+1;
               xa(i,j)=x(count3,j);
               ya(i,j)=y(count3,j);
               za(i,j)=z(count3,j);
               for k=1:10
                   if j==1
                       MSF1f(k,i)=MSF1i(k,count3);
                   else
                       MSF2f(k,i)=MSF2i(k,count3);
                   end
               end
           end
        end
    end
       
    
    [ix,jx]=find(~xa);
    [iy,jy]=find(~ya);
    [iz,jz]=find(~za);
    
   for i=1:length(jx)
       xa(ix(i),2^abs(jx(i)-2))=0;
       ya(iy(i),2^abs(jy(i)-2))=0;
       za(iz(i),2^abs(jz(i)-2))=0;
   end
   
    dixa=xa(:,1)-xa(:,2);
    diya=ya(:,1)-ya(:,2);
    diza=za(:,1)-za(:,2);   
               
    rd=sqrt((sum(dixa.*dixa)+sum(diya.*diya)+sum(diza.*diza))/count3);
    for i=1:length(dixa)
        ird(i)=sqrt(dixa(i)^2+diya(i)^2+diza(i)^2);
    end   
       
     dvs=ird;
end
    
ndvs=dvs/trapz(dvs);
%ndv=dv;

%plot(ndvs,'LineWidth',3)
%hold on
%line([0 981],[mean(ndv(1:981)) mean(ndv(1:981))])
%xline(960)
%xline(960+948)
resnum=length(ndvs);
end
