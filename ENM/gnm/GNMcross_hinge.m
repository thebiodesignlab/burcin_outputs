%% This code is used to find hinge residues from GNM cross-correlations.
clear all; close all; clc;
fname1='1OMP';
chain='A';
modeset=[1];
[cross,resnum]=GNMcross(fname1,modeset,chain)

imagesc(cross,[-1,1]);
set(gca,'YDir','normal','Fontsize',24);
set(gcf,'Color',[1 1 1])
caxis([-1 1]);
colorbar;
axis square
colormap jet

xlabel('Residue number','Fontsize',30)
ylabel('Residue number','Fontsize',30)


if length(modeset)==1
coloring=cross(1,:);

crossfile=strcat('cross',regexprep(int2str(modeset), '  ', '_'),'.txt');
fid = fopen(crossfile,'wt');
fprintf(fid,'%f\n',coloring)
fclose(fid);

count=0;
init=1;
f=0;
for i=2:resnum
    if round(coloring(i)*1000)/1000~=round(coloring(i-1)*1000)/1000
        f=f+1;
        count=count+1;
        F(f,:)=[init i-1 i-init round(coloring(i-1)*1000)/1000];
        hinge(count)=i;
        hinge(count+1)=i-1;
        count=count+1;
        init=i;
    end
end
F(f+1,:)=[init,i,i-init+1, round(coloring(i)*1000)/1000];
lF=F(:,3)>15;
count=0;
ch=0;
for i=1:length(lF)
    j=lF(i);
    if j==0
        G(count,:)=[F(i-1,1) F(i+1,2) F(i+1,2)-F(i-1,1) F(i-1,4)];
        ch=i+1;
    else
        
        if i~=ch
            count=count+1;
            G(count,:)=F(i,:);
        end
    end
end
        
    
hinges=unique(hinge);
hingetxt=strcat('hinges',int2str(modeset),'.txt');
fid=fopen(hingetxt,'wt');
fprintf(fid,' %d',hinges)
fclose(fid);
end

