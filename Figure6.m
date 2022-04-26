% Bang et al (2021) Neurocomputational mechanisms of confidence in self and
% other
%
% Reproduces Figure 6
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
dirDataFMRI_timeSeries= [repoBase,fs,'Data',fs,'fMRI',fs,'ROI_TimeSeries'];
dirDataFMRI_canonHRF= [repoBase,fs,'Data',fs,'fMRI',fs,'ROI_CanonHRF'];
dirFigures= [repoBase,fs,'Figures'];

% Add paths
addpath('Functions');

% ROIs
roi_v= {'V5','LIP'};
cWindow= 'dec'; % time window (dec: decision; gam: gamble)

% Subjects
sbj_v= [22:27 29:32 34:42 44:45]; % subject numbers
k= 0; for j= sbj_v; k= k+1; sbj_name{k}= ['mriData_sub_',num2str(j)]; end % subject names

%% -----------------------------------------------------------------------
%% CANONICAL HRF

% Load exclusions
load([dirDataFMRI_canonHRF,fs,'exclusions_',cWindow,'.mat']);

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
    load([dirDataFMRI_canonHRF,fs,cWindow,'_',roi_v{i_roi},'_s',num2str(sbj_v(k)),'.mat']);
    roi_HRF = timeSeries;
        
    % select non-NaN trials
    nonans = output.exclusion{sbj_v(k)}==0&(data.incl==1);
    
    % GLM
    % Full
    idx = nonans==1;
    roi_Z= zscore(roi_HRF(idx));
    soc= (data.self(idx)==0)-.5;
    coh= zscore(data.theta(idx));
    x= [soc; coh; soc.*coh]';
    y= roi_Z';
    b= glmfit(x,y);
    beta_HRF{i_roi}.full(k,:) = b(2:end);
    % Self
    idx = nonans==1 & data.self==1;
    roi_Z= zscore(roi_HRF(idx));
    coh= zscore(data.theta(idx));
    x= coh';
    y= roi_Z';
    b= glmfit(x,y);
    beta_HRF{i_roi}.self(k,:) = b(2:end);
    % Other
    idx = nonans==1 & data.self==0;
    roi_Z= zscore(roi_HRF(idx));
    coh= zscore(data.theta(idx));
    x= coh';
    y= roi_Z';
    b= glmfit(x,y);
    beta_HRF{i_roi}.other(k,:) = b(2:end);
        
end

end

%% -----------------------------------------------------------------------
%% STATS
% 
% [a b c d]= ttest(beta_HRF{1}.full(:,1));
% MT_soc= b
% [a b c d]= ttest(beta_HRF{1}.full(:,2));
% MT_coh= b
% [a b c d]= ttest(beta_HRF{1}.full(:,3));
% MT_int= b
% [a b c d]= ttest(beta_HRF{1}.self);
% MT_self= b
% [a b c d]= ttest(beta_HRF{1}.other);
% MT_other= b
% 
% [a b c d]= ttest(beta_HRF{1}.full(:,1));
% LIP_soc= b
% [a b c d]= ttest(beta_HRF{1}.full(:,2));
% LIP_coh= b
% [a b c d]= ttest(beta_HRF{2}.full(:,3));
% LIP_int= b
% [a b c d]= ttest(beta_HRF{2}.self);
% LIP_self= b
% [a b c d]= ttest(beta_HRF{2}.other);
% LIP_other= b

%% -----------------------------------------------------------------------
%% VISUALISATION

%% General specifications
scol= [30 144 255]./255;
ocol= [255 165 0]./255;
ccol= [0 1 1];
gcol= [0 1 0];
mcol= [1 0 1];
kcol= [0 0 0];
axisFS= 60;
% labelFS= 44;
alphaz= .5;
lw= 8;
jitter= 1;
ms= 100;

%% FIGURE
% loop through ROIs
for i_roi= 1:length(roi_v);
    % data
    soc= beta_HRF{i_roi}.full(:,1);
    coh= beta_HRF{i_roi}.full(:,2);
    int= beta_HRF{i_roi}.full(:,3);   
    % create figure
    figz=figure('color',[1 1 1]);
    % plot data
    steplotmulticoloffalpha(soc,ccol,lw,1,alphaz);
    steplotmulticoloffalpha(coh,mcol,lw,2,alphaz);
    steplotmulticoloffalpha(int,gcol,lw,3,alphaz);
    plot([0 4],[0 0],'k-','LineWidth',lw);
    % overlay subjects
    % social
    y= soc;
    x= violaPoints(1,soc,jitter);
    scatter(x,y,ms,kcol,'filled','MarkerFaceAlpha',.5);
    % coherence
    y= coh;
    x= violaPoints(2,coh,jitter);
    scatter(x,y,ms,kcol,'filled','MarkerFaceAlpha',.5);
    % interaction
    y= int;
    x= violaPoints(3,int,jitter);
    scatter(x,y,ms,kcol,'filled','MarkerFaceAlpha',.5);
    % tidy up
    ylim([-1 .6]);
    set(gca,'YTick',-.8:.4:.4);
    xlim([0 4]);
    h=gca;
    h.XAxis.Visible = 'off';
    set(gca,'FontSize',axisFS,'LineWidth',lw);
    axis square;
    print('-djpeg','-r300',['Figures',fs,'Figure-6A-HRF-',roi_v{i_roi}]);
end

%% FIGURE
% loop through ROIs
for i_roi= 1:length(roi_v);
    % data
    self= beta_HRF{i_roi}.self(:,1);
    other= beta_HRF{i_roi}.other(:,1);
    % create figure
    figz=figure('color',[1 1 1]);
    % plot data
    steplotmulticoloffalpha(self,scol,lw,1,alphaz);
    steplotmulticoloffalpha(other,ocol,lw,2,alphaz);
    plot([0 3],[0 0],'k-','LineWidth',8);
    % overlay subjects
    % self
    y= self;
    x= violaPoints(1,self,jitter);
    scatter(x,y,ms,kcol,'filled','MarkerFaceAlpha',.5);
    % other
    y= other;
    x= violaPoints(2,other,jitter);
    scatter(x,y,ms,kcol,'filled','MarkerFaceAlpha',.5);
    % tidy up
    ylim([-.4 .4]);
    set(gca,'YTick',-.3:.3:.3);
    xlim([0 3]);
    h=gca;
    h.XAxis.Visible = 'off';
    set(gca,'FontSize',axisFS,'LineWidth',lw);
    axis square;
    print('-djpeg','-r300',['Figures',fs,'Figure-6B-HRF-',roi_v{i_roi}]);
end


%% -----------------------------------------------------------------------
%% TIME COURSE

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
    load([dirDataFMRI_timeSeries,fs,cWindow,'_',roi_v{i_roi},'_s',num2str(sbj_v(k)),'.mat']);
    roi_ts = timeSeries;
    
    % select non-NaN trials
    nonans = ~isnan(roi_ts(:,end))'&(data.incl==1);
    
    % UP-SAMPLED GLM
    % Full
    idx = nonans==1;
    roi_Zts = zscore(roi_ts(idx,:));
    soc= (data.self(idx)==0)-.5;
    coh= zscore(data.theta(idx));
    t= 0;
    for j= 1:size(roi_ts,2)
        t= t+1;
        x= [soc; coh; soc.*coh]';
        y= roi_Zts(:,j);
        beta= glmfit(x,y,'normal');
        beta_ts{i_roi}.soc(k,t) = beta(end-2);
        beta_ts{i_roi}.coh(k,t) = beta(end-1);
        beta_ts{i_roi}.int(k,t) = beta(end-0);
    end
    % Self
    idx = nonans==1 & data.self==1;
    roi_Zts = zscore(roi_ts(idx,:));
    coh= zscore(data.theta(idx));
    t= 0;
    for j= 1:size(roi_ts,2)
        t= t+1;
        x= coh';
        y= roi_Zts(:,j);
        beta= glmfit(x,y,'normal');
        beta_ts{i_roi}.self(k,t) = beta(end);
    end
    % Other
    idx = nonans==1 & data.self==0;
    roi_Zts = zscore(roi_ts(idx,:));
    coh= zscore(data.theta(idx));
    t= 0;
    for j= 1:size(roi_ts,2)
        t= t+1;
        x= coh';
        y= roi_Zts(:,j);
        beta= glmfit(x,y,'normal');
        beta_ts{i_roi}.other(k,t) = beta(end);
    end
    
end

end

%% -----------------------------------------------------------------------
%% VISUALISATION

%% General specifications
scol= [30 144 255]./255;
ocol= [255 165 0]./255;
ccol= [0 1 1];
gcol= [0 1 0];
mcol= [1 0 1];
axisFS= 28;
labelFS= 40;
lw= 4;
max_t = 99; % index for max time
srate = .144; % sampling rate in seconds

%% FIGURE
% loop through ROIs
for i_roi= 1:length(roi_v);  
    % data
    soc= beta_ts{i_roi}.soc;
    coh= beta_ts{i_roi}.coh;
    int= beta_ts{i_roi}.int;
    % create figure
    figz=figure('color',[1 1 1]);
    % add reference lines
    plot([1/srate 1/srate],[-1 +1],'k-','LineWidth',lw/2); hold on
    plot([2.5/srate 2.5/srate],[-1 +1],'k-','LineWidth',lw/2); hold on
    plot([0 max_t],[0 0],'k-','LineWidth',lw); hold on
    % plot beta time series
    fillsteplotcol(soc,lw,'-',ccol); hold on
    fillsteplotcol(coh,lw,'-',mcol); hold on
    fillsteplotcol(int,lw,'-',gcol); hold on
    % tidy up
    ylim([-.54 .54]);
    set(gca,'YTick',-.4:.2:.4);
    xlim([1 84]);
    set(gca,'XTick',[1 14 28 42 56 70 84])
    set(gca,'XTickLabel',{'0','2','4','6','8','10','12'})
    box('off')
    set(gca,'FontSize',axisFS,'LineWidth',lw);
    xlabel('time [seconds]','FontSize',labelFS,'FontWeight','normal');
    ylabel('beta [arb. units]','FontSize',labelFS,'FontWeight','normal');
    print('-djpeg','-r300',['Figures',fs,'Figure-6A-time-',roi_v{i_roi}]);    
end

%% FIGURE
% loop through ROIs
for i_roi= 1:length(roi_v);
    % data
    self= beta_ts{i_roi}.self;
    other= beta_ts{i_roi}.other;    
    % create figure
    figz=figure('color',[1 1 1]);
    % add reference lines
    plot([1/srate 1/srate],[-1 +1],'k-','LineWidth',lw/2); hold on
    plot([2.5/srate 2.5/srate],[-1 +1],'k-','LineWidth',lw/2); hold on
    plot([0 max_t],[0 0],'k-','LineWidth',lw); hold on
    % plot beta time series with significance overlaid
    fillsteplotcol(self,lw,'-',scol); hold on
    fillsteplotcol(other,lw,'-',ocol); hold on
    % tidy up
    ylim([-.28 .28]);
    set(gca,'YTick',-.2:.1:.2);
    xlim([1 84]);
    set(gca,'XTick',[1 14 28 42 56 70 84])
    set(gca,'XTickLabel',{'0','2','4','6','8','10','12'})
    box('off')
    set(gca,'FontSize',axisFS,'LineWidth',lw);
    xlabel('time [seconds]','FontSize',labelFS,'FontWeight','normal');
    ylabel('beta [arb. units]','FontSize',labelFS,'FontWeight','normal');
    print('-djpeg','-r300',['Figures',fs,'Figure-6B-time-',roi_v{i_roi}]);  
end