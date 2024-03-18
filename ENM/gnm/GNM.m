function [MSF,resnum]=GNM(fname1,mode_beg,mode_max,chain)

rcut_gnm=10;
% gamma constant
ga=1;
try
prot=pdbread(fname1);
catch
prot=getpdb(fname1);  
end

atomnum=size(prot.Model.Atom,2);
count=0;
for i=1:atomnum
	if contains(prot.Model.Atom(i).AtomName,'CA')==1 && contains(prot.Model.Atom(i).chainID,chain)==1
        if isempty(prot.Model(1).Atom(i).altLoc) || ...
                        strcmpi(prot.Model(1).Atom(i).altLoc,'A')
	count=count+1;
	x(count)=prot.Model.Atom(i).X;
	y(count)=prot.Model.Atom(i).Y;
	z(count)=prot.Model.Atom(i).Z;
        end
	end
end

resnum=length(x);
if mode_beg<0 && mode_max<0
    mode_beg=resnum+mode_beg;
    mode_max=resnum+mode_max;
end


for j=1:resnum
	for k=1:resnum
         distx = x(j)-x(k);
         disty = y(j)-y(k);
         distz = z(j)-z(k);

         r=sqrt(distx^2+disty^2+distz^2);               
                % Kirchhoff
                if (r <= rcut_gnm && j~=k && r > 0.0001)
                    cont(j,k)=-1;
                else
                    cont(j,k)=0;                                
                end
                
	end
end

        % detailed balance for connectivity matrix
        diagonal=sum(cont(:,:));
        
        for j=1:resnum
            for i=1:resnum
                if i == j
                    cont(i,j)=-1*diagonal(i); % Kirchhoff
                end
            end
        end
    
        [U(:,:),S(:,:),V(:,:)]=svd(cont(:,:));
        w=diag(S(:,:));        
        
%         
%         % GNM 10 slow eigenvectors and slow modes
%         for j=resnum-1:-1:resnum-mode_max
%             for i=1:resnum
%                 slow_vectors_gnm(i,j)=V(i,j);
%                 slow_modes_gnm(i,j)=V(i,j)*V(i,j);
%             end
%         end
        
	
% gnmCfile = fopen('GNM_MSF','a');
invcont=zeros(resnum,resnum);
for kk=resnum-mode_max:resnum-mode_beg
    for j=1:resnum
        for i=1:resnum
            invcont(i,j)=invcont(i,j)+U(i,kk)*U(j,kk)/w(kk);
        end
    end
end

MSF=diag(invcont);
end


