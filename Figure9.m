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
    load([dirDataFMRI_canonHRF,fs,cWindow,'_',roi_v{i_roi},'_s',num2str(sbj_v(k)),'.mat']);
    roi_HRF = timeSeries;
            
    % exclusions
    nonansM= sim.lik<1; % model
    nonansE= data.incl==1; % behaviour
    betaEXL= output.exclusion{sbj_v(k)}==0; % brain
    nonansN= betaEXL(nonansE); % brain
    
    % extraction
    fnames=fieldnames(data); 
    for i_f= 1:length(fnames); evalc(['data.',fnames{i_f},'=data.',fnames{i_f},'(nonansE);']); end
    fnames=fieldnames(sim); 
    for i_f= 1:length(fnames); evalc(['sim.',fnames{i_f},'=sim.',fnames{i_f},'(nonansM);']); end
    roi_HRF= roi_HRF(nonansE);
    
    % GLM
    % Full
    idx = nonansN==1 & data.self==0;
    roi_Z= zscore(roi_HRF(idx));
    acc= data.acc(idx);
    con= sim.ppCo(idx);
    spe= zscore(acc-con);
    x= spe';
    y= roi_Z';
    b= glmfit(x,y);
    beta_HRF{i_roi}.full(k,:) = b(2:end);
    
end

end

% -----------------------------------------------------------------------
% STATS

[a b c d]= ttest(beta_HRF{1}.full(:,1));

[a b c d]= ttest(beta_HRF{2}.full(:,1));

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
    spe= beta_HRF{i_roi}.full(:,1);
    % create figure
    figz=figure('color',[1 1 1]);
    % plot data
    steplotmulticoloffalpha(spe,kcol,lw,1,alphaz);
    plot([0 2],[0 0],'k-','LineWidth',8);
    % overlay subjects
    % spe
    y= spe;
    x= violaPoints(1,spe,jitter);
    scatter(x,y,ms,kcol,'filled','MarkerFaceAlpha',.5);
    % tidy up
    ylim([-.6 .6]);
    set(gca,'YTick',-.4:.4:.4);
    xlim([0 2]);
    h=gca;
    h.XAxis.Visible = 'off';
    set(gca,'FontSize',axisFS,'LineWidth',lw);
    axis square;
    print('-djpeg','-r300',['Figures',fs,'Figure-9-HRF-',roi_v{i_roi}]);
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
    
    % select non-NaN trials
    nonans = ~isnan(roi_ts(:,end))'&(data.incl==1);
    
    % UP-SAMPLED GLM
    % Full
    idx = nonans==1 & data.self==0;
    roi_Zts = zscore(roi_ts(idx,:));
    acc = data.acc(idx);
    con = sim.ppCo(idx);
    pe = zscore(acc-con);
    t= 0;
    for j= 1:size(roi_ts,2)
        t= t+1;
        x= pe';
        y= roi_Zts(:,j);
        beta= glmfit(x,y,'normal');
        beta_ts{i_roi}.pe(k,t) = beta(end);
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
    pe= beta_ts{i_roi}.pe;
    % create figure
    figz=figure('color',[1 1 1]);
    % add reference lines
    plot([2/srate 2/srate],[-1 +1],'k-','LineWidth',lw/2); hold on
    plot([5/srate 5/srate],[-1 +1],'k-','LineWidth',lw/2); hold on
    plot([0 max_t],[0 0],'k-','LineWidth',lw); hold on
    % plot beta time series
    fillsteplotcol(pe,lw,'-',kcol); hold on
    % tidy up
    ylim([-.24 .24]);
    set(gca,'YTick',-.4:.1:.4);
    xlim([1 84]);
    set(gca,'XTick',[1 14 28 42 56 70 84])
    set(gca,'XTickLabel',{'0','2','4','6','8','10','12'})
    box('off')
    set(gca,'FontSize',axisFS,'LineWidth',lw);
    xlabel('time [seconds]','FontSize',labelFS,'FontWeight','normal');
    ylabel('beta [arb. units]','FontSize',labelFS,'FontWeight','normal');
    print('-djpeg','-r300',['Figures',fs,'Figure-9-time-',roi_v{i_roi}]);    
end