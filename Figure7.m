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
    idx = nonansN==1;
    roi_Z= zscore(roi_HRF(idx));
    soc= (data.self(idx)==0)-.5;
    coh= zscore(data.theta(idx));
    conS= sim.ppCs(idx);
    conS(data.self(idx)==0)= 0;
    conO= sim.ppCo(idx);
    conO(data.self(idx)==1)= 0;
    con= conS+conO;
    con= zscore(con);
    skl= zscore(data.kskill(idx));
    orthogonal= spm_orth([coh; con]'); % requires SPM
    con= zscore(orthogonal(:,end))';
    x= [soc; con; soc.*con]';
    y= roi_Z';
    b= glmfit(x,y);
    beta_HRF{i_roi}.full(k,:) = b(2:end);
    % Self
    idx = nonansN==1 & data.self==1;
    roi_Z= zscore(roi_HRF(idx));
    coh= zscore(data.theta(idx));
    con= zscore(sim.ppCs(idx));
    orthogonal= spm_orth([coh; con]'); % requires SPM
    con= zscore(orthogonal(:,end))';
    x= con';
    y= roi_Z';
    b= glmfit(x,y);
    beta_HRF{i_roi}.self(k,:) = b(2:end);
    % Other
    idx = nonansN==1 & data.self==0;
    roi_Z= zscore(roi_HRF(idx));
    coh= zscore(data.theta(idx));
    con= zscore(sim.ppCo(idx));
    orthogonal= spm_orth([coh; con]'); % requires SPM
    con= zscore(orthogonal(:,end))';
    x= con';
    y= roi_Z';
    b= glmfit(x,y);
    beta_HRF{i_roi}.other(k,:) = b(2:end);
 
end

end

%% -----------------------------------------------------------------------
%% STATS
% 
% [a b c d]= ttest(beta_HRF{1}.full(:,1));
% TPJ_soc= b
% [a b c d]= ttest(beta_HRF{1}.full(:,2));
% TPJ_con= b
% [a b c d]= ttest(beta_HRF{1}.full(:,3));
% TPJ_int= b
% [a b c d]= ttest(beta_HRF{1}.self);
% TPJ_self= b
% [a b c d]= ttest(beta_HRF{1}.other);
% TPJ_other= b
% 
% [a b c d]= ttest(beta_HRF{2}.full(:,1));
% dmPFC_soc= b
% [a b c d]= ttest(beta_HRF{2}.full(:,2));
% dmPFC_con= b
% [a b c d]= ttest(beta_HRF{2}.full(:,3));
% dmPFC_int= b
% [a b c d]= ttest(beta_HRF{2}.self);
% dmPFC_self= b
% [a b c d]= ttest(beta_HRF{2}.other);
% dmPFC_other= b


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

%% General specifications
kcol= [0 0 0];
axisFS= 60;
% labelFS= 44;
alphaz= .5;
lw= 8;

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
    plot([0 4],[0 0],'k-','LineWidth',8);
    % tidy up
    ylim([-.48 .48]);
    set(gca,'YTick',-.4:.4:.4);
    xlim([0 4]);
    h=gca;
    h.XAxis.Visible = 'off';
    set(gca,'FontSize',axisFS,'LineWidth',lw);
    axis square;
    print('-djpeg','-r300',['Figures',fs,'Figure-7A-HRF-',roi_v{i_roi}]);
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
    % tidy up
    ylim([-.16 .16]);
    set(gca,'YTick',-.1:.1:.1);
    xlim([0 3]);
    h=gca;
    h.XAxis.Visible = 'off';
    set(gca,'FontSize',axisFS,'LineWidth',lw);
    axis square;
    print('-djpeg','-r300',['Figures',fs,'Figure-7B-HRF-',roi_v{i_roi}]);
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
    idx = nonans==1;
    roi_Zts = zscore(roi_ts(idx,:));
    soc= (data.self(idx)==0)-.5;
    coh= zscore(data.theta(idx));
    conS= sim.ppCs(idx);
    conS(data.self(idx)==0)= 0;
    conO= sim.ppCo(idx);
    conO(data.self(idx)==1)= 0;
    con= conS+conO;
    con= zscore(con);
    skl= zscore(data.kskill(idx));
    orthogonal= spm_orth([coh; con]'); % requires SPM
    con= zscore(orthogonal(:,end))';
    t= 0;
    for j= 1:size(roi_ts,2)
        t= t+1;
        x= [soc; con; soc.*con]';
        y= roi_Zts(:,j);
        beta= glmfit(x,y,'normal');
        beta_ts{i_roi}.soc(k,t) = beta(end-2);
        beta_ts{i_roi}.con(k,t) = beta(end-1);
        beta_ts{i_roi}.int(k,t) = beta(end-0);
    end
    % Self
    idx = nonans==1 & data.self==1;
    roi_Zts = zscore(roi_ts(idx,:));
    coh= zscore(data.theta(idx));
    con= zscore(sim.ppCs(idx));
    orthogonal= spm_orth([coh; con]'); % requires SPM
    con= zscore(orthogonal(:,end))';
    t= 0;
    for j= 1:size(roi_ts,2)
        t= t+1;
        x= con';
        y= roi_Zts(:,j);
        beta= glmfit(x,y,'normal');
        beta_ts{i_roi}.self(k,t) = beta(end);
    end
    % Other
    idx = nonans==1 & data.self==0;
    roi_Zts = zscore(roi_ts(idx,:));
    coh= zscore(data.theta(idx));
    con= zscore(sim.ppCo(idx));
    orthogonal= spm_orth([coh; con]'); % requires SPM
    con= zscore(orthogonal(:,end))';
    t= 0;
    for j= 1:size(roi_ts,2)
        t= t+1;
        x= con';
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
    con= beta_ts{i_roi}.con;
    int= beta_ts{i_roi}.int;
    % create figure
    figz=figure('color',[1 1 1]);
    % add reference lines
    plot([1/srate 1/srate],[-1 +1],'k-','LineWidth',lw/2); hold on
    plot([2.5/srate 2.5/srate],[-1 +1],'k-','LineWidth',lw/2); hold on
    plot([0 max_t],[0 0],'k-','LineWidth',lw); hold on
    % plot beta time series
    fillsteplotcol(soc,lw,'-',ccol); hold on
    fillsteplotcol(con,lw,'-',mcol); hold on
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
    ylabel('beta [a.u.]','FontSize',labelFS,'FontWeight','normal');
    print('-djpeg','-r300',['Figures',fs,'Figure-7A-time-',roi_v{i_roi}]);    
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
    ylabel('beta [a.u.]','FontSize',labelFS,'FontWeight','normal');
    print('-djpeg','-r300',['Figures',fs,'Figure-7B-time-',roi_v{i_roi}]);  
end