clear all; clc;close all;
addpath '/Users/acarb/Dropbox (The Francis Crick)/DeBenedictisE/burcin/crick_biodesign/postdoc/trymatlab'
%%Input Data
fname1='1OMP';
fname2='1ANF';
chain1='A';
chain2='A';
%% Calculate Dif Vec
[resnum,ndvs]=difvecPDB(fname1,chain1,fname2,chain2);
mode_max=10;
mode_beg=1;
% coarseDV=smoothdata(ndvs,'movmean',4); 
%lapDV=del2(coarseDV);
% ndvs01=ndvs;
% TFn=islocalmin(ndvs,'FlatSelection','all');
% ndvs01(TFn)=1;
% ndvs01(~TFn)=0;
windowSize = 1; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
ndvsf=filter(b,a,ndvs);
% figure(1)
% plot(ndvs)
% hold on
% plot(ndvsf)
%% Calculate MSF and overlaps
count=0
for kk=resnum-mode_max:resnum-mode_beg
    count=count+1;
    [MSF,resnum]=GNM(fname1,1,count,chain1);
    MSFf=filter(b,a,MSF);
%     MSF01=MSF;
%     TFm=islocalmin(MSF,'FlatSelection','all');
%     MSF01(TFm)=1;
%     MSF01(~TFn)=0;
    ovM(1,count)=dot(MSFf,ndvsf)/(norm(MSFf)*norm(ndvsf)); % regular overlap
    %coarseMSF=smoothdata(MSF,'movmean',4); 
    %ovM(2,count)=dot(coarseMSF,coarseDV)/(norm(coarseMSF)*norm(coarseDV));
   % ovM(3,count)=dot(MSF01,ndvs01)/(norm(MSF01)*norm(ndvs01));
%     c=corrcoef(MSF,ndvs);
%     ovM(3,count)=c(2,1);
    %lapMSF=del2(coarseMSF);
    %ovM(3,count)=dot(lapMSF,lapDV)/(norm(lapMSF)*norm(lapDV));
     figure(1)
     subplot(2,5,count)
     plot(MSFf,'LineWidth',3)
     hold on
     plot(ndvsf,'LineWidth',3)
%     figure(2)
%     subplot(2,5,count)
%     plot(lapMSF,'LineWidth',3)
%     hold on
%     plot(lapDV,'LineWidth',3)
end
% figure(3)
% imagesc(mode_beg:mode_max,1:3,ovM)
% set(gca,'YDir','normal','Fontsize',24);
% colormap('jet')
% colorbar


