clear all; clc; close all;
addpath '/Users/acarb/Dropbox (The Francis Crick)/DeBenedictisE/burcin/crick_biodesign/postdoc/trymatlab'
%[MSF(1,:),resnum]=GNM('1EU8',1,10,'A');
for i=1:10
[MSF(i,:),resnum]=GNM('1EU8',i,i,'A');

fileID = fopen('BFAC_enrich_r0r3.txt','r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);

x0=MSF(i,:);
y0=A(2:length(MSF)+1);
nz=find(y0);
x=x0(nz);
y=y0(nz);
p = polyfit(x, y, 1);
px = [min(x) max(x)];
py = polyval(p, px);
subplot(2,5,i)
%subplot(1,1,i)
scatter(x,y,130,'filled','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','k')
hold on
plot(px, py, ':','LineWidth', 3);
title(['Mode ',num2str(i)])
ylabel('Enrichment Score')
xlabel('Fluctuation')
set(gca,'FontSize',30)
axis square
grid on
xlim([min(x) max(x)])
end
