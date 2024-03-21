%% ANM code by Burcin Acar. Given a PDB file named fname1, this code
% calculates mode_max number of slow ANM modes starting from the 1st.
% It outputs MSFs of each individual mode.
function [nanmcross]=anmcross(fname1,mode_max,chain)

try
    prot=getpdb(fname1);
catch
    prot=pdbread(fname1);  
end
atomnum=size(prot.Model.Atom,2);
count=0;

for i=1:atomnum
	if contains(prot.Model.Atom(i).AtomName,'CA')==1&& contains(prot.Model.Atom(i).chainID,chain)==1 
        if isempty(prot.Model(1).Atom(i).altLoc) || ...
                        strcmpi(prot.Model(1).Atom(i).altLoc,'A')
            count=count+1;
            x(count)=prot.Model.Atom(i).X;
            y(count)=prot.Model.Atom(i).Y;
            z(count)=prot.Model.Atom(i).Z;
        end
    end
end


%  general
resnum=size(x,2); 
anmcross=zeros(resnum,resnum);
nanmcross=anmcross;
rcut_anm= 13;

rcut2= rcut_anm*rcut_anm; %???

resnum3=resnum*3;
hessian=zeros(resnum3,resnum3);
% gamma constant
ga=1;

        
        for j=1:resnum
            for k=1:resnum
                distx = x(1,j)-x(1,k);   
                disty = y(1,j)-y(1,k);   
                distz = z(1,j)-z(1,k);
                
                r2(j,k)=distx^2+disty^2+distz^2;
                
                r=sqrt(distx^2+disty^2+distz^2);
                
                
                % Hessian
                if (r <= rcut_anm && k~=j)
                    % Creation of Hii
                    % diagonals
                    hessian(3*j-2,3*j-2)=hessian(3*j-2,3*j-2)+ga*distx*distx/r2(j,k); % x, 3j-2
                    hessian(3*j-1,3*j-1)=hessian(3*j-1,3*j-1)+ga*disty*disty/r2(j,k); % y, 3j-1
                    hessian(3*j,3*j)=hessian(3*j,3*j)+ga*distz*distz/r2(j,k);         % z, 3j
                    
                    % off-diagonals
                    hessian(3*j-2,3*j-1)=hessian(3*j-2,3*j-1)+ga*distx*disty/r2(j,k); % xy
                    hessian(3*j-2,3*j)=hessian(3*j-2,3*j)+ga*distx*distz/r2(j,k);     % xz
                    
                    hessian(3*j-1,3*j-2)=hessian(3*j-1,3*j-2)+ga*disty*distx/r2(j,k); % yx
                    hessian(3*j-1,3*j)=hessian(3*j-1,3*j)+ga*disty*distz/r2(j,k);     % yz
                    
                    hessian(3*j,3*j-2)=hessian(3*j,3*j-2)+ga*distz*distx/r2(j,k); % zx
                    hessian(3*j,3*j-1)=hessian(3*j,3*j-1)+ga*distz*disty/r2(j,k); % zy
                                        
                    % Creation of Hij
                    % diagonals
                    hessian(3*j-2,3*k-2)=-1*ga*distx*distx/r2(j,k); % x, 3j-2
                    hessian(3*j-1,3*k-1)=-1*ga*disty*disty/r2(j,k); % y, 3j-1
                    hessian(3*j,3*k)=    -1*ga*distz*distz/r2(j,k); % z, 3j
                    
                    % off-diagonals
                    hessian(3*j-2,3*k-1)=-1*ga*distx*disty/r2(j,k); % xy
                    hessian(3*j-2,3*k)=  -1*ga*distx*distz/r2(j,k); % xz
                    
                    hessian(3*j-1,3*k-2)=-1*ga*disty*distx/r2(j,k); % yx
                    hessian(3*j-1,3*k)=  -1*ga*disty*distz/r2(j,k); % yz
                    
                    hessian(3*j,3*k-2)=  -1*ga*distz*distx/r2(j,k); % zx
                    hessian(3*j,3*k-1)=  -1*ga*distz*disty/r2(j,k); % zy

                end
            end
        end
		
        
        [Ua(:,:),Sa(:,:),Va(:,:)]=svd(hessian(:,:));
        wa=diag(Sa(:,:));
        % to export as output
        w_snsh_a(:,1)=diag(Sa(:,:));    
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % ANM 20 slow eigenvectors and slow modes, last column corresponds to first mode
        for j=resnum3-6:-1:1%resnum3-(mode_max+5) % resnum3-25 for 20 modes
            for i=1:resnum3
                slow_vectors_anm(i,j)=Va(i,j);

            end
        end

        
        % ANM
        for m=1:mode_max%20
            for n=1:resnum
                Vx(n,m)=slow_vectors_anm(3*n-2,resnum3-5-m);
                Vy(n,m)=slow_vectors_anm(3*n-1,resnum3-5-m);
                Vz(n,m)=slow_vectors_anm(3*n  ,resnum3-5-m);    
            end
        end

        % ANM first column corresponds to first mode
        for m=1:mode_max%20
            
             for j=1:resnum
           
                slow_modes_anm(j,m)=(Vx(j,m)*Vx(j,m)+Vy(j,m)*Vy(j,m)+Vz(j,m)*Vz(j,m))/wa(resnum3-5-m);
                slow_modes_not_weighted_anm(j,m)=sqrt((Vx(j,m)*Vx(j,m)+Vy(j,m)*Vy(j,m)+Vz(j,m)*Vz(j,m)));
                for i=1:resnum
                    anmcross(i,j)=anmcross(i,j)+(Vx(j,m)*Vx(i,m)+Vy(j,m)*Vy(i,m)+Vz(j,m)*Vz(i,m))/wa(resnum3-5-m);
                    
                end
            end
        end
        for i=1:resnum
            for j=1:resnum
                nanmcross(i,j)=anmcross(i,j)/sqrt(anmcross(i,i)*anmcross(j,j));
            end
        end
        imagesc(nanmcross,[-1,1]);
        set(gca,'YDir','normal','Fontsize',24);
        set(gcf,'Color',[1 1 1])
        caxis([-1 1]);
        colorbar;
        axis square
        colormap jet

        xlabel('Residue number','Fontsize',30)
        ylabel('Residue number','Fontsize',30)
end

