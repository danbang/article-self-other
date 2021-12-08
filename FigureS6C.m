% Bang et al (2021) Neurocomputational mechanisms of confidence in self and
% other
%
% Reproduces Figure S6B
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
bestModel= 'T_B2T2P1_U';
Nreruns= 4;

% loop through reruns
for i_r= 1:Nreruns;
    load([dirModelFits,fs,'model_',bestModel,'_seed_',num2str(i_r),'_vb.mat']);
    tmp_group_param_MAP(i_r,:)= group_param_MAP;
end

% average group-level parameters
mu_group_param_MAP= mean(tmp_group_param_MAP);

% input for visualisation of utility function
valueL= fliplr(0:.01:1);
valueG= 0:.01:1;
powerL= mu_group_param_MAP(1);
powerG= mu_group_param_MAP(2);
averse= mu_group_param_MAP(3);

% loop through values
j= 0;
% loss domain
for i_L= 1:length(valueL);
    j= j+1;
    utility(j)= -averse*(valueL(i_L)^powerL);
    value(j)= -valueL(i_L);
end
% gain domain
for i_G= 1:length(valueG);
    j= j+1;
    utility(j)= (valueG(i_G)^powerG);
    value(j)= valueG(i_G);
end

%% -----------------------------------------------------------------------
%% VISUALISATION

% create figure
xlims= [-1 1];
ylims= [-3 3];
figz=figure('color',[1 1 1]);
hold on;
plot([xlims(1) xlims(2)],[-1.5 -1.5],'k-','LineWidth',2);
plot([xlims(1) xlims(2)],[0 0],'k-','LineWidth',4);
plot([xlims(1) xlims(2)],[1.5 1.5],'k-','LineWidth',2);
plot([-.5 -.5],[ylims(1) ylims(2)],'k-','LineWidth',2);
plot([0 0],[ylims(1) ylims(2)],'k-','LineWidth',4);
plot([.5 .5],[ylims(1) ylims(2)],'k-','LineWidth',2);
plot(value,utility,'m-','LineWidth',4);
set(gca,'FontSize',28,'LineWidth',4);
set(gca,'XTick',-1:.5:1,'YTick',[-3:1:3]);
xlabel('value','FontSize',40);
ylabel('utility','FontSize',40);
xlim(xlims);
ylim(ylims);
box on;
axis square;
print('-djpeg','-r300',['Figures',fs,'Figure-S6C']);