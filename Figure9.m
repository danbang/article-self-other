% Bang et al (2021) Neurocomputational mechanisms of confidence in self and
% other
%
% Reproduces Figure 9
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
dirModelFits= [repoBase,fs,'Modelling',fs,'Stan_fits'];
dirFigures= [repoBase,fs,'Figures'];

% Add paths
addpath('Functions');

% Specify winning model for overlaying onto behavioural data
bestModel= 'T_B2T2P1';
Nreruns= 4;

% ROIs
roi_v= {'TPJ','dmPFC'};
cWindow= 'gam'; % time window (dec: decision; gam: gamble)

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
    
    % load model data
    clear sim tmp;
    % loop through reruns
    for i_r= 1:Nreruns;
        load([dirModelFits,fs,'model_',bestModel,'_seed_',num2str(i_r),'_vb.mat']);
        tmp= subject_predict;
        fnames=fieldnames(tmp); 
        for i_f= 1:length(fnames); evalc(['sim.',fnames{i_f},'(i_r,:)=tmp.',fnames{i_f},'(k,:);']); end
    end
    % average over seeds
    fnames=fieldnames(tmp); 
    for i_f= 1:length(fnames); evalc(['sim.',fnames{i_f},'=mean(sim.',fnames{i_f},');']); end
    
    % load neural data
    load([dirDataFMRI,fs,cWindow,'_',roi_v{i_roi},'_s',num2str(sbj_v(k)),'.mat']);
    roi_ts = timeSeries;
    
    % exclusions (M: model; E: behavioural; N: neural)
    nonansM= sim.lik<1;
    nonansE= data.incl==1;
    roi_ts = roi_ts(nonansE,:);
    nonansN= ~isnan(roi_ts(:,end))';
    
    % extraction
    fnames=fieldnames(data); 
    for i_f= 1:length(fnames); evalc(['data.',fnames{i_f},'=data.',fnames{i_f},'(nonansE);']); end
    fnames=fieldnames(sim); 
    for i_f= 1:length(fnames); evalc(['sim.',fnames{i_f},'=sim.',fnames{i_f},'(nonansM);']); end

    % UP-SAMPLED GLM
    nonansN = nonansN & data.self==0;
    roi_Zts = zscore(roi_ts(nonansN,:));
    acc     = data.acc(nonansN);
    con     = sim.ppCs(nonansN);
    pe      = zscore(acc-con);
    t= 0;
    for j= 1:size(roi_ts,2)
        t= t+1;
        x= [pe]';
        y= roi_Zts(:,j);
        beta= glmfit(x,y,'normal','constant','on');
        beta_ts{i_roi}.pe(k,t) = beta(end);
    end
    
end

end

%% Load permutation statistics
load('Permutation/Figure9.mat');

%% -----------------------------------------------------------------------
%% VISUALISATION

%% General specifications
kcol= [0 0 0];
axisFS= 28;
labelFS= 40;
lw= 4;
max_t = 99; % index for max time
srate = .144; % sampling rate in seconds

%% FIGURE: ENCODING COHERENCE
% loop through ROIs
for i_roi= 1:length(roi_v);
    
    % Permutation test
    % data
    pe= beta_ts{i_roi}.pe;
    % t-tests
    [H,P,CI,STATS]= ttest(pe);
    tstat_pe= STATS.tstat;
    % 95% bounds
    tstat_permutation_LB_pe= quantile(permutation_beta_ts{i_roi}.pe,.025);
    tstat_permutation_UB_pe= quantile(permutation_beta_ts{i_roi}.pe,.975);
    % Compare to bounds
    tstat_significant_pe= (tstat_pe<tstat_permutation_LB_pe)|(tstat_pe>tstat_permutation_UB_pe);
    
    % create figure
    figz=figure('color',[1 1 1]);
    % add reference lines
    plot([2/srate 2/srate],[-1 +1],'k-','LineWidth',lw/2); hold on
    plot([5/srate 5/srate],[-1 +1],'k-','LineWidth',lw/2); hold on
    plot([0 max_t],[0 0],'k-','LineWidth',lw); hold on
    % plot beta time series with significance overlaid
    fillsteplotcol(pe,lw,'-',kcol); hold on
    p = tstat_significant_pe; for i= 1:length(p); if p(i); plot(i,-.17,'s','Color',kcol,'MarkerFaceColor',kcol); hold on; end; end;
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
    print('-djpeg','-r300',['Figures',fs,'Figure-9-',roi_v{i_roi}]);
end