% Bang et al (2021) Neurocomputational mechanisms of confidence in self and
% other
%
% Reproduces Figure 7B
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
    idx= data.self(nonansN)==0;
    roi= roi_ts(nonansN,:);
    roi= zscore(roi(idx,:));
    skl= zscore(data.kskill(idx));
    coh= zscore(data.theta(idx));
    con= zscore(sim.ppCs(idx));
    orthogonal= spm_orth([coh; con]'); % requires SPM
    con= zscore(orthogonal(:,end))';
    t= 0;
    for j= 1:size(roi,2)
        t= t+1;
        x= [con]';
        y= roi(:,j);
        beta= glmfit(x,y,'normal','constant','on');
        beta_ts{i_roi}.con(k,t) = beta(end);
    end
    
end

end


%% General specifications
gcol= [0 1 0];
mcol= [1 0 1];
axisFS= 28;
labelFS= 40;
lw= 4;
max_t = 99; % index for max time
srate = .144; % sampling rate in seconds

%% FIGURE: ENCODING CONFIDENCE
% create figure
figz=figure('color',[1 1 1]);
% add reference lines
plot([1/srate 1/srate],[-1 +1],'k-','LineWidth',lw/2); hold on
plot([2.5/srate 2.5/srate],[-1 +1],'k-','LineWidth',lw/2); hold on
plot([0 max_t],[0 0],'k-','LineWidth',lw); hold on
% plot beta time series with significance overlaid
% TPJ
trace= beta_ts{1}.con;
fillsteplotcol(trace,lw,'-',mcol); hold on
p = ttest(trace); for i= 1:length(p); if p(i); plot(i,-.12,'s','Color',mcol,'MarkerFaceColor',mcol); hold on; end; end;
% dmPFC
trace= beta_ts{2}.con;
fillsteplotcol(trace,lw,'-',gcol); hold on
p = ttest(trace); for i= 1:length(p); if p(i); plot(i,-.13,'s','Color',gcol,'MarkerFaceColor',gcol); hold on; end; end;
% tidy up
ylim([-.15 .15]);
set(gca,'YTick',-.2:.1:.2);
xlim([1 max_t]);
set(gca,'XTick',1:14:max_t-1)
set(gca,'XTickLabel',{'0','2','4','6','8','10','12'})
box('off')
set(gca,'FontSize',axisFS,'LineWidth',lw);
xlabel('time [seconds]','FontSize',labelFS,'FontWeight','normal');
ylabel('beta [a.u.]','FontSize',labelFS,'FontWeight','normal');
print('-djpeg','-r300',['Figures',fs,'Figure7B']);