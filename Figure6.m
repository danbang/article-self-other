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
dirDataFMRI= [repoBase,fs,'Data',fs,'fMRI',fs,'ROI_TimeSeries'];
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
    nonans = ~isnan(roi_ts(:,end))'&(data.incl==1);
    
    % UP-SAMPLED GLM
    % Self
    idx = nonans==1 & data.self==1;
    roi_Zts = zscore(roi_ts(idx,:));
    type= data.self(idx)-.5;
    coh= zscore(data.theta(idx));
    rt1= zscore(data.rt1(idx));
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
    type= data.self(idx)-.5;
    coh= zscore(data.theta(idx));
    rt1= zscore(data.rt1(idx));
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
axisFS= 28;
labelFS= 40;
lw= 4;
max_t = 99; % index for max time
srate = .144; % sampling rate in seconds

%% FIGURE: ENCODING COHERENCE
% loop through ROIs
for i_roi= 1:length(roi_v);
    % create figure
    figz=figure('color',[1 1 1]);
    % add reference lines
    plot([1/srate 1/srate],[-1 +1],'k-','LineWidth',lw/2); hold on
    plot([2.5/srate 2.5/srate],[-1 +1],'k-','LineWidth',lw/2); hold on
    plot([0 max_t],[0 0],'k-','LineWidth',lw); hold on
    % plot beta time series with significance overlaid
    self= beta_ts{i_roi}.self;
    other= beta_ts{i_roi}.other;
    fillsteplotcol(self,lw,'-',scol); hold on
    p = ttest(self); for i= 1:length(p); if p(i); plot(i,-.17,'s','Color',scol,'MarkerFaceColor',scol); hold on; end; end;
    fillsteplotcol(other,lw,'-',ocol); hold on
    p = ttest(other); for i= 1:length(p); if p(i); plot(i,-.18,'s','Color',ocol,'MarkerFaceColor',ocol); hold on; end; end;
    % tidy up
    ylim([-.21 .21]);
    set(gca,'YTick',-.2:.1:.2);
    xlim([1 84]);
    set(gca,'XTick',[1 14 28 42 56 70 84])
    set(gca,'XTickLabel',{'0','2','4','6','8','10','12'})
    box('off')
    set(gca,'FontSize',axisFS,'LineWidth',lw);
    xlabel('time [seconds]','FontSize',labelFS,'FontWeight','normal');
    ylabel('beta [a.u.]','FontSize',labelFS,'FontWeight','normal');
    print('-djpeg','-r300',['Figures',fs,'Figure6_',roi_v{i_roi}]);
end