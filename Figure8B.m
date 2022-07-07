% Bang et al (2021) Neurocomputational mechanisms of confidence in self and
% other
%
% Reproduces Figure 8B
%
% Dan Bang danbang.db@gmail.com 2021

%% -----------------------------------------------------------------------
%% PREPARATION

% Fresh memory
clear; close all;

% Paths [change 'repoBase' according to local setup]
fs= filesep;
repoBase= [getDropbox(1),fs,'Ego',fs,'Matlab',fs,'ucl',fs,'self_other',fs,'Repository'];
dirDataBehaviour= [repoBase,fs,'Data',fs,'Behaviour',fs,'Scan',fs,'Task'];
dirDataFMRI= [repoBase,fs,'Data',fs,'fMRI',fs,'ROI_PPI'];
dirFigures= [repoBase,fs,'Figures'];

% Add paths
addpath('Functions');

% ROIs
seed= 'LIP';
roi_v= {'TPJ','dmPFC'};
cWindow= 'dec'; % time window (dec: decision; gam: gamble)

%% -----------------------------------------------------------------------
%% ANALAYSIS

% load PPI seed
load([dirDataFMRI,fs,'ppi','_',cWindow,'_',seed,'_SvsO.mat']);

% log ROI PPI contrast estimates
for i_roi = 1:length(roi_v)
    
    % load ROI PPI contrast estimates
    ppi(:,i_roi)= -ppiConEst{i_roi}; % negative: SvsO -> OvsS

end

%% -----------------------------------------------------------------------
%% VISUALISATION

%% General specifications
colz= [1 0 1;
       0 1 0];
dcol= [0 0 0];
alphaz= .6;
lw= 4;
ms= 24;
axisFS= 28;
labelFS= 40;
jitter=.5;

%% FIGURE: PPI CONTRAST estimates
% create figure
figure('color',[1 1 1]);
% plot group mean and SEM
for i_roi= 1:size(ppi,2);
        steplotmulticoloffalpha(ppi(:,i_roi),colz(i_roi,:),lw,i_roi,alphaz);
end
% overlay subjects
for i=1:size(ppi,2);
    y= ppi(:,i);
    x= violaPoints(i,ppi(:,i),jitter);
    scatter(x,y,ms,dcol,'filled','MarkerFaceAlpha',.4);
end
% add reference line
plot([0 size(ppi,2)+1],[0 0],'k-','LineWidth',lw);
% tidy up
set(gca,'Xtick',[1:length(roi_v)],'XTickLabel',roi_v,'XTickLabelRotation',0);
set(gca,'FontSize',axisFS,'LineWidth',lw);
xlabel('ROI','FontWeight','normal','FontSize',labelFS);
ylabel('PPI connectivity','FontWeight','normal','FontSize',labelFS);
ylim([-0.6 +0.6]); 
set(gca,'YTick',-.6:.2:.6);
xlim([0 length(roi_v)+1]);
box('off')
print('-djpeg','-r300',['Figures',fs,'Figure-8B']);

% write .csv source data file
csvheaders= {'panel_b_TPJ','panel_b_dmPFC'};
csvdata= [ppi(:,1) ppi(:,2)];
csvwrite_with_headers('SourceData/Figure8',csvdata,csvheaders);