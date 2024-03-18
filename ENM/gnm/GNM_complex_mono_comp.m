%% This code compares the modes derived from the complex with the
% allosteric modes of the monomer.
%% Initialize
clear all; clc; close all;
set(0,'DefaultFigureWindowStyle','docked')
monomodeset=[1 7];
%monof='1OMP'; % apo structure 
monof='1ANF';
complexf1='3OSQ';
complexf2='3OSR';

[dummy,MSF1,resnum]=GNM(monof,monomodeset(1),'A');
[dummy,MSF2,resnum]=GNM(monof,monomodeset(2),'A');
clear dummy
% merge 173 & 427 for 3OSQ
% merge 311 & 564 for 3OSR 

[MSFc1,dummy,resnum2]=GNM(complexf1,1:10,'A');
[MSFc2,dummy,resnum3]=GNM(complexf2,1:10,'A');
clear dummy

for kk=1:10
    MSFc1m(kk,:)=MSFc1(kk,[1:173-1 422-6-1-3:resnum2]);
    MSFc2m(kk,:)=MSFc2(kk,[5:311+4 561-17-3+4:resnum3]);
    ov11(kk)=dot(MSFc1m(kk,:),MSF1(2:369))/(norm(MSFc1m(kk,:))*norm(MSF1(2:369)));
    ov13(kk)=dot(MSFc1m(kk,:),MSF2(2:369))/(norm(MSFc1m(kk,:))*norm(MSF2(2:369)));
    ov21(kk)=dot(MSFc2m(kk,:),MSF1(1:370))/(norm(MSFc2m(kk,:))*norm(MSF1(1:370)));
    ov23(kk)=dot(MSFc2m(kk,:),MSF2(1:370))/(norm(MSFc2m(kk,:))*norm(MSF2(1:370)));
end


[i11 j11]=max(ov11); % Complex 1 resembling monomodeset(1)
[i13 j13]=max(ov13); % Complex 1 resembling monomodeset(2)
[i21 j21]=max(ov21); % Complex 2 resembling monomodeset(1)
[i23 j23]=max(ov23); % Complex 2 resembling monomodeset(2)



[dummy,MSFavg,resnum2]=GNM(complexf1,sort([j11 j13]),'A'); %complex 1 avg
[dummy,MSFavgm,resnum]=GNM(monof,monomodeset,'A'); %monomer avg
MSFc1mavg=MSFavg([1:173-1 422-6-1-3:resnum2]); 
[dummy,MSFavg2,resnum3]=GNM(complexf2,sort([j21 j23]),'A'); %complex 2 avg
MSFc2mavg=MSFavg([5:311+4 561-17-3+4:resnum3]); 
 %% averages 
figure(40)
plot(2:369,MSFavgm(2:369),'LineWidth',3)
hold on
plot(2:369,MSFc1mavg./trapz(MSFc1mavg),'LineWidth',3)
set(gca,'FontSize',24)
xlabel('Residue Number')
grid on
axis square
ylabel('Mode Fluctuation')
legend('monomer','GFP-fused')
title('3OSQ')

figure(41)
plot(1:370,MSFavgm(1:370),'LineWidth',3)
hold on
plot(1:370,MSFc2mavg./trapz(MSFc2mavg),'LineWidth',3)
set(gca,'FontSize',24)
xlabel('Residue Number')
grid on
axis square
ylabel('Mode Fluctuation')
legend('monomer','GFP-fused')
title('3OSR')
%% individual modes
figure(1)
plot(2:369,MSF1(2:369),'LineWidth',3)
hold on
plot(2:369,MSFc1m(j11,:)./trapz(MSFc1m(j11,:)),'LineWidth',3)
set(gca,'FontSize',24)
xlabel('Residue Number')
grid on
axis square
ylabel('Mode Fluctuation')
legend('monomer1','GFP-fused1')
title('GNM Mode',j11)

figure(2)
plot(2:369,MSF2(2:369),'LineWidth',3)
hold on
plot(2:369,MSFc1m(j13,:)./trapz(MSFc1m(j13,:)),'LineWidth',3)
set(gca,'FontSize',24)
xlabel('Residue Number')
grid on
axis square
ylabel('Mode Fluctuation')
legend('monomer2','GFP-fused1')
title('GNM Mode',j13)

figure(3)
plot(1:370,MSF1,'LineWidth',3)
hold on
plot(1:370,MSFc2m(j21,:)./trapz(MSFc2m(j21,:)),'LineWidth',3)
set(gca,'FontSize',24)
xlabel('Residue Number')
grid on
axis square
ylabel('Mode Fluctuation')
legend('monomer1','GFP-fused2')
title('GNM Mode',j21)

figure(4)
plot(1:370,MSF2,'LineWidth',3)
hold on
plot(1:370,MSFc2m(j23,:)./trapz(MSFc2m(j23,:)),'LineWidth',3)
set(gca,'FontSize',24)
xlabel('Residue Number')
grid on
axis square
ylabel('Mode Fluctuation')
legend('monomer2','GFP-fused2')
title('GNM Mode',j23)



