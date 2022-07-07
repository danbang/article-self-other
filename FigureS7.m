% Bang et al (2021) Neurocomputational mechanisms of confidence in self and
% other
%
% Reproduces Figure S4
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
dirDataFMRI= [repoBase,fs,'Data',fs,'fMRI',fs,'ROI_CanonHRF'];
dirFigures= [repoBase,fs,'Figures'];

% Add paths
addpath('Functions');

% ROIs
roi_v= {'V5','LIP','TPJ','dmPFC'};
cWindow= 'dec'; % time window (dec: decision; gam: gamble)

% Load exclusions
load([dirDataFMRI,fs,'exclusions_',cWindow,'.mat']);

% Subjects
sbj_v= [22:27 29:32 34:42 44:45]; % subject numbers
k= 0; for j= sbj_v; k= k+1; sbj_name{k}= ['mriData_sub_',num2str(j)]; end % subject names

%% -----------------------------------------------------------------------
%% ANALAYSIS

% loop through ROIs
for i_roi = 1:length(roi_v)

% initialise subject index
k = 0;

% loop through subjects
for i_sbj = 1:length(sbj_v)
        
    % update subject index
    k=k+1;
    
    % load behavioural data
    load([dirDataBehaviour,fs,sbj_name{k},'.mat']);
    concatenate;
    fnames=fieldnames(tmp); 
    for i_f= 1:length(fnames); evalc(['data.',fnames{i_f},'=tmp.',fnames{i_f},';']); end
    
    % load neural data
    load([dirDataFMRI,fs,cWindow,'_',roi_v{i_roi},'_s',num2str(sbj_v(k)),'.mat']);
    roi_ts = timeSeries;
    
    % select non-NaN trials
    nonans = output.exclusion{sbj_v(k)}==0&(data.incl==1);
        
    % further selection because we are looking at trial pairs
    incl= ones(1,120);
    for t= 1:120;
        if data.trial(t)<2;
            incl(t)= 0;
        end
    end
    for t=1:120
        if nonans(t)==0;
            incl(t)= 0;
        end
    end
    for t=2:120
        if nonans(t-1)==0;
            incl(t)= 0;
        end
    end
    
    % construct variables
    indx= incl==1;
    slfA= data.self;
    slfC= slfA;
    slfP= [0 slfA(1:end-1)];
    slfC= slfC(indx);
    slfP= slfP(indx);
    roiC= zscore(roi_ts(indx));
    
    % perform analysis
    for r= 1:2;
        for c= 1:2;
            meanz{i_roi}(r,c,k)= mean(roiC( slfC==r-1 & slfP==c-1 ));
        end      
    end
    
    % regression
    x= [zscore(slfC); zscore(slfP);zscore(slfC).*zscore(slfP)]';
    y= roiC';
    b= glmfit(x,y);
    betaz{i_roi}(k,:) = b(2:end);
    
end

end


%% -----------------------------------------------------------------------
%% VISUALISATION

%% General specifications
scol= [30 144 255]./255;
ocol= [255 165 0]./255;
dcol= [0 0 0];
alphaz= .6;
lw= 4;
ms= 20;
axisFS= 24;
labelFS= 36;
jitter=.5;
my_labels= {'MT+','LIP','TPJ','dmPFC'};

%% FIGURE: REPETITION SUPPRESSION
% loop through ROIs
for i_roi= 1:length(roi_v)
    % create figure
    figz=figure('color',[1 1 1]);
    % select data
    OO= squeeze(meanz{i_roi}(1,1,:));
    SO= squeeze(meanz{i_roi}(1,2,:));
    OS= squeeze(meanz{i_roi}(2,1,:));
    SS= squeeze(meanz{i_roi}(2,2,:));
    tmpz= [OO SO OS SS];
    % plot group data
    steplotmulticoloffalpha(OO,ocol,lw,1,.4);
    steplotmulticoloffalpha(SO,ocol,lw,2,.8);
    steplotmulticoloffalpha(OS,scol,lw,3,.4);
    steplotmulticoloffalpha(SS,scol,lw,4,.8);
    % overlay subject data
    for i=1:4;
        y= tmpz(:,i);
        x= violaPoints(i,tmpz(:,i),jitter);
        scatter(x,y,ms,dcol,'filled','MarkerFaceAlpha',.4);
    end
    % add reference line
    plot([0 5],[0 0],'k-','LineWidth',lw)
    % tidy up
    set(gca,'XTick',[1:4],'XTickLabel',{'O -> O','S -> O','O -> S','S -> S'});
    set(gca,'YTick',[-.6:.3:.6]);
    ylim([-.70 .70]); 
    xlim([0 5]);
    box('off')
    set(gca,'FontSize',axisFS,'LineWidth',lw);
    xlabel('condition','FontSize',labelFS,'FontWeight','normal');
    ylabel('activity [z-score]','FontSize',labelFS,'FontWeight','normal');
    title(my_labels{i_roi},'Fontsize',labelFS,'FontWeight','normal')
    print('-djpeg','-r300',['Figures',fs,'Figure-S7-',roi_v{i_roi}]);
end

% write .csv source data file
csvheaders= {'V5_O2O','V5_S2O','V5_02S','V5_S2S', ...
             'LIP_O2O','LIP_S2O','LIP_02S','LIP_S2S', ...
             'TPJ_O2O','TPJ_S2O','TPJ_02S','TPJ_S2S', ...
             'dmPFC_O2O','dmPFC_S2O','dmPFC_02S','dmPFC_S2S'};
csvdata= [squeeze(meanz{1}(1,1,:)) squeeze(meanz{1}(1,2,:)) squeeze(meanz{1}(2,1,:)) squeeze(meanz{1}(2,2,:)), ...
          squeeze(meanz{2}(1,1,:)) squeeze(meanz{2}(1,2,:)) squeeze(meanz{2}(2,1,:)) squeeze(meanz{2}(2,2,:)), ...
          squeeze(meanz{3}(1,1,:)) squeeze(meanz{3}(1,2,:)) squeeze(meanz{3}(2,1,:)) squeeze(meanz{3}(2,2,:)), ...
          squeeze(meanz{4}(1,1,:)) squeeze(meanz{4}(1,2,:)) squeeze(meanz{4}(2,1,:)) squeeze(meanz{4}(2,2,:))];
csvwrite_with_headers('SourceData/FigureS7',csvdata,csvheaders);