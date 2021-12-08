% Bang et al (2021) Neurocomputational mechanisms of confidence in self and
% other
%
% Reproduces Figure 7
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
    idx= nonansN;
    soc= zscore(data.self(idx))*-1;
    roi= roi_ts(idx,:);
    roi= zscore(roi(idx,:),1,1);
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
    for j= 1:size(roi,2)
        t= t+1;
        x= [soc; con; soc.*con]';
        y= roi(:,j);
        beta= glmfit(x,y,'normal','constant','on');
        beta_ts{i_roi}.soc(k,t) = beta(2);
        beta_ts{i_roi}.con(k,t) = beta(3);
        beta_ts{i_roi}.int(k,t) = beta(4);
    end
    
end

end

%% Load permutation statistics
load('Permutation/Figure7.mat');

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

%% FIGURE: GLM
for i_roi= 1:2;
    
    % Permutation test
    % data
    soc= beta_ts{i_roi}.soc;
    con= beta_ts{i_roi}.con;
    int= beta_ts{i_roi}.int;
    % t-tests
    [H,P,CI,STATS]= ttest(soc);
    tstat_soc= STATS.tstat;
    [H,P,CI,STATS]= ttest(con);
    tstat_con= STATS.tstat;
    [H,P,CI,STATS]= ttest(int);
    tstat_int= STATS.tstat;
    % 95% bounds
    tstat_permutation_LB_soc= quantile(permutation_beta_ts{i_roi}.slf,.025);
    tstat_permutation_UB_soc= quantile(permutation_beta_ts{i_roi}.slf,.975);
    tstat_permutation_LB_con= quantile(permutation_beta_ts{i_roi}.con,.025);
    tstat_permutation_UB_con= quantile(permutation_beta_ts{i_roi}.con,.975);
    tstat_permutation_LB_int= quantile(permutation_beta_ts{i_roi}.int,.025);
    tstat_permutation_UB_int= quantile(permutation_beta_ts{i_roi}.int,.975);
    % Compare to bounds
    tstat_significant_soc= (tstat_soc<tstat_permutation_LB_soc)|(tstat_soc>tstat_permutation_UB_soc);
    tstat_significant_con= (tstat_con<tstat_permutation_LB_con)|(tstat_con>tstat_permutation_UB_con);
    tstat_significant_int= (tstat_int<tstat_permutation_LB_int)|(tstat_int>tstat_permutation_UB_int);
    
    % create figure
    figz=figure('color',[1 1 1]);
    % add reference lines
    plot([1/srate 1/srate],[-1 +1],'k-','LineWidth',lw/2); hold on
    plot([2.5/srate 2.5/srate],[-1 +1],'k-','LineWidth',lw/2); hold on
    plot([0 max_t],[0 0],'k-','LineWidth',lw); hold on
     % plot beta time series with significance overlaid
    fillsteplotcol(soc,lw,'-',ccol); hold on
    p = tstat_significant_soc; for i= 1:length(p); if p(i); plot(i,-.16,'s','Color',ccol,'MarkerFaceColor',ccol); hold on; end; end;
    fillsteplotcol(con,lw,'-',mcol); hold on
    p = tstat_significant_con; for i= 1:length(p); if p(i); plot(i,-.18,'s','Color',mcol,'MarkerFaceColor',mcol); hold on; end; end;
    fillsteplotcol(int,lw,'-',gcol); hold on
    p = tstat_significant_int; for i= 1:length(p); if p(i); plot(i,-.20,'s','Color',gcol,'MarkerFaceColor',gcol); hold on; end; end;
    % tidy up
    ylim([-.28 .28]);
    set(gca,'YTick',-.3:.1:.3);
    xlim([1 84]);
    set(gca,'XTick',[1 14 28 42 56 70 84])
    set(gca,'XTickLabel',{'0','2','4','6','8','10','12'})
    box('off')
    set(gca,'FontSize',axisFS,'LineWidth',lw);
    xlabel('time [seconds]','FontSize',labelFS,'FontWeight','normal');
    ylabel('beta [a.u.]','FontSize',labelFS,'FontWeight','normal');
    % save figure
    print('-djpeg','-r300',['Figures',fs,'Figure-7A-',roi_v{i_roi}]);
    
    % log overlap in significance for visualisation (ignore late time points)
    p_overlap2{i_roi}= (tstat_significant_soc(1:70)+tstat_significant_con(1:70)+tstat_significant_int(1:70))==2;
    p_overlap3{i_roi}= (tstat_significant_soc(1:70)+tstat_significant_con(1:70)+tstat_significant_int(1:70))==3;
    
end

%% FIGURE: VISUALISATION
for i_roi= 1:2;
    
    % data
    soc= beta_ts{i_roi}.soc;
    con= beta_ts{i_roi}.con;
    int= beta_ts{i_roi}.int;
    
    % ROI-specific
    if i_roi== 1;
        current_timepoints= find(p_overlap3{i_roi}==1);
    elseif i_roi== 2;
        current_timepoints= find(p_overlap2{i_roi}==1);
    end
    
    % average across subjects and then timepoints
    muBeta_soc= mean(mean(beta_ts{i_roi}.soc(:,current_timepoints)));
    muBeta_con= mean(mean(beta_ts{i_roi}.con(:,current_timepoints)));
    muBeta_int= mean(mean(beta_ts{i_roi}.int(:,current_timepoints)));
    
    % vectors for calculating predicted activity
    v_soc= [-1 1];
    v_con= -3:.1:3;

    % calculate predicted activity
    for i_soc= 1:length(v_soc);
        for i_con= 1:length(v_con);
            activity(i_soc,i_con)= muBeta_soc*v_soc(i_soc) + muBeta_con*v_con(i_con) + muBeta_int*v_soc(i_soc)*v_con(i_con);
        end
    end
    
    % create figure
    figz=figure('color',[1 1 1]);
    hold on;
    % plot reference line
    plot([-3.5 3.5],[0 0],'k-','LineWidth',lw);
    % plot predictions
    plot(v_con,activity(1,:),'-','color',scol,'LineWidth',lw);
    plot(v_con,activity(2,:),'-','color',ocol,'LineWidth',lw); 
    % tidy up
    ylim([-.45 .45]);
    set(gca,'YTick',-.4:.2:.4);
    xlim([-3.5 3.5]);
    set(gca,'XTick',-3:1:3)
    box('off')
    set(gca,'FontSize',axisFS,'LineWidth',lw);
    xlabel('confidence','FontSize',labelFS,'FontWeight','normal');
    ylabel('estimated activity','FontSize',labelFS,'FontWeight','normal');
    % save figure
    print('-djpeg','-r300',['Figures',fs,'Figure-7B-',roi_v{i_roi}]);
    
end