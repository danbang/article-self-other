% Bang et al (2021) Neurocomputational mechanisms of confidence in self and
% other
%
% Reproduces Figure S6B
%
% Dan Bang danbang.db@gmail.com 2021

%% -----------------------------------------------------------------------
%% PREPARATION

% Fresh memory
clear; %close all;

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
ToM_Zero= 'T_B2T2P1';
ToM_Util= 'T_B2T2P1_U';
Nreruns= 4;

% Subjects
sbj_v= [22:27 29:32 34:42 44:45]; % subject numbers
k= 0; for j= sbj_v; k= k+1; sbj_name{k}= ['mriData_sub_',num2str(j)]; end % subject names

%% -----------------------------------------------------------------------
%% ANALAYSIS

k= 0;

% loop through subjects
for i_sbj = 1:length(sbj_v)
        
    % update subject index
    k=k+1;
    
    % load behavioural data
    load([dirDataBehaviour,fs,sbj_name{k},'.mat']);
    concatenate;
    fnames=fieldnames(tmp); 
    for i_f= 1:length(fnames); evalc(['data.',fnames{i_f},'=tmp.',fnames{i_f},';']); end
    
    % load ToM_Zero data
    clear sim tmp;
    % loop through reruns
    for i_r= 1:Nreruns;
        load([dirModelFits,fs,'model_',ToM_Zero,'_seed_',num2str(i_r),'_vb.mat']);
        tmp= subject_predict;
        fnames=fieldnames(tmp); 
        for i_f= 1:length(fnames); evalc(['sim.',fnames{i_f},'(i_r,:)=tmp.',fnames{i_f},'(k,:);']); end
    end
    % average over seeds
    fnames=fieldnames(tmp); 
    for i_f= 1:length(fnames); evalc(['sim.',fnames{i_f},'=mean(sim.',fnames{i_f},');']); end 
    % relabel
    sim_ToM_Zero= sim;
    
    % load ToM_util data
    clear sim tmp;
    % loop through reruns
    for i_r= 1:Nreruns;
        load([dirModelFits,fs,'model_',ToM_Util,'_seed_',num2str(i_r),'_vb.mat']);
        tmp= subject_predict;
        fnames=fieldnames(tmp); 
        for i_f= 1:length(fnames); evalc(['sim.',fnames{i_f},'(i_r,:)=tmp.',fnames{i_f},'(k,:);']); end
    end
    % average over seeds
    fnames=fieldnames(tmp); 
    for i_f= 1:length(fnames); evalc(['sim.',fnames{i_f},'=mean(sim.',fnames{i_f},');']); end
    % relabel
    sim_ToM_Util= sim;
    
    % exclusions (M: model; E: behavioural; N: neural)
    nonansM= sim.lik<1;
    nonansE= data.incl==1;
    
    % re-extraction
    fnames=fieldnames(data); 
    for i_f= 1:length(fnames); evalc(['data.',fnames{i_f},'=data.',fnames{i_f},'(nonansE);']); end
    fnames=fieldnames(sim_ToM_Zero); 
    for i_f= 1:length(fnames); evalc(['sim_ToM_Zero.',fnames{i_f},'=sim_ToM_Zero.',fnames{i_f},'(nonansM);']); end
    fnames=fieldnames(sim_ToM_Util); 
    for i_f= 1:length(fnames); evalc(['sim_ToM_Util.',fnames{i_f},'=sim_ToM_Util.',fnames{i_f},'(nonansM);']); end

    % correlations
    % variables
    con_Zero_self= sim_ToM_Zero.ppCs(data.self==1);
    con_Zero_other= sim_ToM_Zero.ppCo(data.self==0);
    con_Util_self= sim_ToM_Util.ppCs(data.self==1);
    con_Util_other= sim_ToM_Util.ppCo(data.self==0);
    con_Zero_both= [con_Zero_self con_Zero_other];
    con_Util_both= [con_Util_self con_Util_other];
    % glm
    output= fitlm(con_Zero_self,con_Util_self);
    rsquare_con_self(i_sbj,1)= output.Rsquared.Ordinary;
    output= fitlm(con_Zero_other,con_Util_other);
    rsquare_con_other(i_sbj,1)= output.Rsquared.Ordinary;
    output= fitlm(con_Zero_both,con_Util_both);
    rsquare_con_both(i_sbj,1)= output.Rsquared.Ordinary;
    % pearson's
    [r p]= corrcoef(con_Zero_self,con_Util_self);
    r_con_self(i_sbj,1)= r(2);
    [r p]= corrcoef(con_Zero_other,con_Util_other);
    r_con_other(i_sbj,1)= r(2);
    output= corrcoef(con_Zero_both,con_Util_both);
    r_con_both(i_sbj,1)= r(2);
    
end

%% -----------------------------------------------------------------------
%% VISUALISATION

%% General specifications
scol= [30 144 255]./255;
ocol= [255 165 0]./255;
scol2= [30 144 255]./295;
ocol2= [255 165 0]./295;
dcol= [0 0 0];
alphaz= .6;
lw= 4;
ms= 200;
axisFS= 24;
labelFS= 36;
jitter=.5;

%% FIGURE: GAMBLE BY COHERENCE
% Create figure
figure('color',[1 1 1]);
hold on;
for i_x= 1:21;
    plot([i_x i_x],[0 1],'k-');
end
y_v= .1:.1:1;
for i_y= 1:length(y_v);
    plot([0 22],[y_v(i_y) y_v(i_y)],'k-');
end
scatter((1:length(sbj_v))-.1,r_con_self,ms,scol,'filled','MarkerFaceAlpha',.7);
scatter((1:length(sbj_v))+.1,r_con_other,ms,ocol,'filled','MarkerFaceAlpha',.7);
% plot(rsquare_con_self,'ko','MarkerFaceColor',scol,'MarkerSize',ms)
% plot(rsquare_con_other,'ko','MarkerFaceColor',ocol,'MarkerSize',ms)
% Tidy up
xlim([0 22]);
ylim([0 1]); 
set(gca,'XTick',1:21,'XTickLabel',{'1','2','3','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','21'});
set(gca,'YTick',0:.1:1);
set(gca,'LineWidth',lw,'FontSize',axisFS);
ylabel('\it{r}','FontSize',labelFS);
xlabel('subject','FontSize',labelFS);
box on;
print('-djpeg','-r300',['Figures/Figure-S6B']);