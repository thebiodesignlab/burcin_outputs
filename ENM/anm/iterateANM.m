% This code looks for allostery as coupling between different functional sites, defined by the cross-correlation
% derived from ANM usinf anmcross script.

clear all; clc; close all;

mode_max=10;
chain='A';
%% Sequences to be fused
MBP='KIEEGKLVIWINGDKGYNGLAEVGKKFEKDTGIKVTVEHPDKLEEKFPQVAATGDGPDIIFWAHDRFGGYAQSGLLAEITPDKAFQDKLYPFTWDAVRYNGKLIAYPIAVEALSLIYNKDLLPNPPKTWEEIPALDKELKAKGKSALMFNLQEPYFTWPLIAADGGYAFKYENGKYDIKDVGVDNAGAKAGLTFLVDLIKNKHMNADTDYSIAEAAFNKGETAMTINGPWAWSNIDTSKVNYGVTVLPTFKGQPSKPFVGVLSAGINAASPNKELAKEFLENYLLTDEGLEAVNKDKPLGAVALKSYEEELAKDPRIAATMENAQKGEIMPNIPQMSAFWYAVRTAVINAASGRQTVDEALKDAQTRITK'
cpegfp='SYNVFIMADKQKNGIKANFKIRHNIEDGSVQLAYHYQQNTPIGDGPVLLPDNHYLSVQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYKGGTGGSMVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYIQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFN'
linker1='VAS';
linker2='AST'; 
%MBPN,linker1,cpegfp,linker2,MBPC 

%% Original Res No for Functional sites
F1anf=[14, 15, 62, 65, 66, 111, 153, 155, 12, 63, 154, 230, 340];
F4eul=[42  44  46  60  61  62  63  64  66  68  69  92  94  96 112 121 145 146 ...
 148 150 165 167 203 205 220 222]; %enhanced gfp functional sites
a=F4eul(F4eul<=147);
a=a+92+6;
b=F4eul(F4eul>147);
b=b-147;
Fcp4eul=[a b]+1; %circularly permuted enhanced gfp functional sites
sumcross=zeros(1,length(MBP)+1);
sumposcross=sumcross;
%% Update Res No for Each Fusion
%for i=0:length(MBP) 
for i=169:169 
    F1anf(F1anf>i)=F1anf(F1anf>i)+length(linker1)+length(cpegfp)+length(linker2);
    F1(i+1,:)=F1anf;
    F1anf=[14, 15, 62, 65, 66, 111, 153, 155, 12, 63, 154, 230, 340];
    F2(i+1,:)=Fcp4eul+i+length(linker1);
    %fname1=sprintf('fusion%u/fusion%u/ranked_0.pdb',i,i);
    fname1='../ranked_0.pdb';
    ac{i+1}=anmcross(fname1,mode_max,chain);
  
    
    for k=F1(i+1,:)
        for l=F2(i+1,:)
            sumcross(1,i+1)=sumcross(1,i+1)+ac{i+1}(k,l);
            if ac{i+1}(k,l)>0
                sumposcross(1,i+1)=sumposcross(1,i+1)+ac{i+1}(k,l);
            end
        end
    end
    
    
end
% 
% plot(0:length(MBP),sumcross)
% plot(0:1,sumcross(1:2))
% plot(0:1,sumposcross(1:2))
% plot(0:length(MBP),sumposcross)
% figure(2)
% cross2=ac{2}
% imagesc(cross2,[-1,1]);
% set(gca,'YDir','normal','Fontsize',24);
% set(gcf,'Color',[1 1 1])
% caxis([-1 1]);
% colorbar;
% axis square
% colormap jet

% xlabel('Residue number','Fontsize',30)
% ylabel('Residue number','Fontsize',30)
% save('MBP_GFP_nomin')
