%function [sauc]=cluster_ROC(xlsxfilepath,matfile)
%% This script by Burcin Acar, reads seq enrichment data from excel file
%  and functional coupling predictions from a mat file to calculate
%  ROC curve and clustering of functional data.

%% Reading Input
C = readvars("../41467_2016_BFncomms12266_MOESM717_ESM.xlsx","Sheet","en","Range","E5:E374");
%C = readvars(xlsxfilepath,"Sheet","en","Range","E5:E374");
%load('MBP_GFP_nomin.mat',"MBP","sumposcross"); % .mat file including whole
% cross-correlations.
%load('../small.mat'); % .mat file including only MBP and sumposcross
%load(matfile,"MBP","sumposcross")
matfile='sasa.mat';
load(matfile,"MBP","sasa_output2")
sumposcross=sasa_output2(:,2)';
X(:,1)=1:length(MBP);
X(:,2)=sumposcross';
X(:,3)=C;
%X(:,4)=kmedoids(sumposcross', 2); % use cluster medians method
X(:,4)=kmeans(sumposcross', 2); % use cluster mean method

%% How to Cluster Observed Data
dummy=nan(size(X(:,4)));
dummy(C>0)=1;
dummy(C<0)=2;
dummy(C==0)=2;
X(:,5)=dummy;


%% Total No of Trues and Falses
sumP=sum(dummy(dummy==1));
sumN=sum(dummy(dummy==2))/2;
%sumnan=sum(isnan(dummy));

%% Rank predicted data without NaN elements
B=X(:,2);
[k,l]=sort(B(~isnan(dummy)),'descend');
M=dummy(~isnan(dummy));
d=zeros(size(k));
sump=zeros(length(k),1);
sumn=sump;
tpr=sump;
fpr=sumn;
auc=zeros(length(k)-1,1);
d(k>=0.1)=1;
d(k<0.1)=2;



for i=1:length(k)
    for j=1:i
        sump(i)=sump(i)+sum(d(j)==1 && d(j)==M(l(j)));
        sumn(i)=sumn(i)+sum(d(j)==1 && d(j)~=M(l(j)));
    end
    if i>1
        sump(i)= sump(i)+sump(i-1)
        sumn(i)= sumn(i)+sumn(i-1)
    end
    tpr(i)=sump(i)./sumP; %tp/p
    fpr(i)=sumn(i)./sumN; %fp/n
end
for i=1:length(k)-1
    auc(i)=(fpr(i)-fpr(i+1))*(tpr(i));
end
sauc=sum(auc);

%% ROC Curve
figure(1)
plot(fpr,tpr,'LineWidth',3)
xlim([0 1])
ylim([0 1])
set(gca,'FontSize',30)
xlabel('False Positive Rate')
ylabel('True Positive Rate')
grid on
axis square

% %% Cluster Plot
% color = 'br';
% figure(2)
% cnt=0;
% h=zeros(length(MBP),1);
% for i=1:length(MBP)
%     h(i)=scatter(X(i,1),X(i,3),36,color(X(i,4)),'filled'); 
%     hold on
%     if i>1 && cnt==0 && X(i-1,4)~=X(i,4)
%         if X(i-1,4)==1
%             leg=[h(i-1) h(i)];
%             cnt=cnt+1;
%         else
%             leg=[h(i) h(i-1)];
%             cnt=cnt+1;
%         end
%             
%     end
% end
% set(gca,'FontSize',30)
% xlabel('Residue Number')
% ylabel('Enrichment Score')
% grid on
% axis square
% legend(leg,'Positive Prediction','Negative Prediction')
%end