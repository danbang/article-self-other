% Bang et al (2021) Neurocomputational mechanisms of confidence in self and
% other
%
% Reproduces Figure 8C
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
roi_v= {'LIP','TPJ','dmPFC'};
cWindow= 'dec'; % time window (dec: decision; gam: gamble)

% Subjects
sbj_v= [22:27 29:32 34:42 44:45]; % subject numbers
k= 0; for j= sbj_v; k= k+1; sbj_name{k}= ['mriData_sub_',num2str(j)]; end % subject names

%% -----------------------------------------------------------------------
%% ANALAYSIS

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
    load([dirDataFMRI,fs,cWindow,'_',roi_v{1},'_s',num2str(sbj_v(k)),'.mat']);
    LIP_ts = timeSeries;
    load([dirDataFMRI,fs,cWindow,'_',roi_v{2},'_s',num2str(sbj_v(k)),'.mat']);
    TPJ_ts = timeSeries;
    load([dirDataFMRI,fs,cWindow,'_',roi_v{3},'_s',num2str(sbj_v(k)),'.mat']);
    dmPFC_ts = timeSeries;
    
    % select non-NaN trials
    nonans = ~isnan(LIP_ts(:,end))'&(data.incl==1);
        
    % UP-SAMPLED GLM: OTHER
    % index
    idx= nonans;
    % roi
    LIP= zscore(LIP_ts(idx,:));
    TPJ= zscore(TPJ_ts(idx,:));
    dmPFC= zscore(dmPFC_ts(idx,:));
    % task
    soc= (data.self(idx)==0)-.5;
    coh= zscore(data.theta(idx));
    % psychological variable and covariates
    psy = soc;
    cov = coh;
    t= 0;
    % loop through time
    for j= 1:size(LIP,2)
        t= t+1;
        % original analysis
        % LIP-TPJ
        x= [psy' LIP(:,j) LIP(:,j).*psy'];
        y= TPJ(:,j);
        B= glmfit(x,y);
        ppi{1}.X(i_sbj,t)= B(end);  
        % LIP-dmPFC
        x= [psy' LIP(:,j) LIP(:,j).*psy'];
        y= dmPFC(:,j);
        B= glmfit(x,y);
        ppi{2}.X(i_sbj,t)= B(end);
%         % control for motion coherence
%         % LIP-TPJ
%         x= [cov' LIP(:,j).*cov' psy' LIP(:,j) LIP(:,j).*psy']; % control for motion coherence
%         y= TPJ(:,j);
%         B= glmfit(x,y);
%         ppi{1}.psy(i_sbj,t)= B(end);  
%         % LIP-dmPFC
%         x= [cov' LIP(:,j).*cov' psy' LIP(:,j) LIP(:,j).*psy']; % control for motion coherence
%         y= dmPFC(:,j);
%         B= glmfit(x,y);
%         ppi{2}.psy(i_sbj,t)= B(end);
    end
        
end


%% -----------------------------------------------------------------------
%% VISUALISATION

%% General specifications
gcol= [0 1 0];
mcol= [1 0 1];
axisFS= 28;
labelFS= 40;
lw= 4;
max_t = 99; % index for max time
srate = .144; % sampling rate in seconds

%% FIGURE: PPI
% data
TPJ= ppi{1}.X;
dmPFC= ppi{2}.X;
% create figure
figz=figure('color',[1 1 1]);
% add reference line
plot([1/srate 1/srate],[-1 +1],'k-','LineWidth',lw/2); hold on
plot([2.5/srate 2.5/srate],[-1 +1],'k-','LineWidth',lw/2); hold on
plot([0 max_t],[0 0],'k-','LineWidth',lw); hold on
% plot beta time series with significance overlaid
% TPJ
fillsteplotcol(TPJ,lw,'-',mcol); hold on
% dmPFC
fillsteplotcol(dmPFC,lw,'-',gcol); hold on
% tidy up
ylim([-.25 .25]);
set(gca,'YTick',-.2:.1:.2);
xlim([1 84]);
set(gca,'XTick',[1 14 28 42 56 70 84])
set(gca,'XTickLabel',{'0','2','4','6','8','10','12'})
box('off')
set(gca,'FontSize',axisFS,'LineWidth',lw);
xlabel('time [seconds]','FontSize',labelFS,'FontWeight','normal');
ylabel('beta [arb. units]','FontSize',labelFS,'FontWeight','normal');
print('-djpeg','-r300',['Figures',fs,'Figure-8C']);