% Bang et al (2021) Neurocomputational mechanisms of confidence in self and
% other
%
% Reproduces Figure 7A
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
dirModelFits= [repoBase,fs,'Modelling',fs,'Stan_fits'];
dirFigures= [repoBase,fs,'Figures'];

% Add paths
addpath('Functions');

% Specify example model and subject
exampleModel= 'T_B2T2P1';
exampleSubject= 14;
Nreruns= 4;

% Subjects
sbj_v= [22:27 29:32 34:42 44:45]; % subject numbers
k= 0; for j= sbj_v; k= k+1; sbj_name{k}= ['mriData_sub_',num2str(j)]; end % subject names


%% -----------------------------------------------------------------------
%% ANALAYSIS

% initialise subject index
k = 0;

% loop through subjects
for i_s= 1:size(sbj_v,2)

    % update subject index
    k=k+1;

    % load behavioural data
    load([dirDataBehaviour,fs,sbj_name{k},'.mat']);
    concatenate;
    fnames=fieldnames(tmp); 
    for i_f= 1:length(fnames); evalc(['data.',fnames{i_f},'=tmp.',fnames{i_f},';']); end

    % load model data
    % initialise
    my_lik= [];
    my_ppk= [];
    my_ppx= [];
    my_ppCs= [];
    my_ppCo= [];
    my_ppPgamble= [];
    my_ppd= [];
    my_ppu= [];
    my_nan= [];
    % concatenate reruns
    for i_r= 1:Nreruns
        load([dirModelFits,fs,'model_',exampleModel,'_seed_',num2str(i_r),'_vb.mat']);
        nonanz= subject_predict.lik(i_s,:)<1;
        my_lik= [my_lik; subject_predict.lik(i_s,nonanz)];
        my_ppk= [my_ppk; subject_predict.ppk(i_s,nonanz)];
        my_ppx= [my_ppx; subject_predict.ppx(i_s,nonanz)];
        my_ppCs= [my_ppCs; subject_predict.ppCs(i_s,nonanz)];
        my_ppCo= [my_ppCo; subject_predict.ppCo(i_s,nonanz)];
        my_ppPgamble= [my_ppPgamble; subject_predict.ppPgamble(i_s,nonanz)];
        my_ppd= [my_ppd; subject_predict.ppd(i_s,nonanz)];
        my_ppu= [my_ppu; subject_predict.ppu(i_s,nonanz)];
        my_nan= [my_nan; ~nonanz];
    end
    % average across reruns
    sim{i_s}.self= data.self;
    sim{i_s}.skill= data.kskill;
    sim{i_s}.lik= nanmean(my_lik);
    sim{i_s}.k= nanmean(my_ppk);
    sim{i_s}.x= nanmean(my_ppx);
    sim{i_s}.Cs= nanmean(my_ppCs);
    sim{i_s}.Co= nanmean(my_ppCo);
    sim{i_s}.Pgamble= nanmean(my_ppPgamble);
    sim{i_s}.d= nanmean(my_ppd);
    sim{i_s}.u= nanmean(my_ppu);
    sim{i_s}.nan= my_nan(1,:);

end

%% -----------------------------------------------------------------------
%% VISUALISATION

%% General specifications
s= 14;
lw= 4;
ms= 20;
axisFS= 28;
labelFS= 40;
col1= [.8 .8 .8];
col2= [.5 .5 .5];
col3= [.1 .1 .1];


%% FIGURE: LATENT SENSORY NOISE ESTIMATE
% create figure
figure('color',[1 1 1]);
% plot model data
s= exampleSubject;
plot(sim{s}.k(sim{s}.skill==1&sim{s}.self==0),'-','Color',col1,'LineWidth',lw); hold on;
plot(sim{s}.k(sim{s}.skill==2&sim{s}.self==0),'-','Color',col2,'LineWidth',lw); hold on;
plot(sim{s}.k(sim{s}.skill==3&sim{s}.self==0),'-','Color',col3,'LineWidth',lw); hold on;
% tidy up
set(gca,'LineWidth',lw,'FontSize',axisFS);
xlabel('trial','FontSize',labelFS);
ylabel('others'' sensory noise','FontSize',labelFS);
set(gca,'XTick',linspace(2,20,7));
set(gca,'YTick',[0:.1:.5]);
ylim([0 .55]);
xlim([0 21]);
box off;
print('-djpeg','-r300',['Figures/Figure7A_noise']);

%% FIGURE: LATENT CONFIDENCE
% create figure
figure('color',[1 1 1]);
% plot model data
plot(sim{s}.Co(sim{s}.skill==1&sim{s}.self==0),'-','Color',col1,'LineWidth',lw); hold on;
plot(sim{s}.Co(sim{s}.skill==2&sim{s}.self==0),'-','Color',col2,'LineWidth',lw); hold on;
plot(sim{s}.Co(sim{s}.skill==3&sim{s}.self==0),'-','Color',col3,'LineWidth',4); hold on;
% tidy up
set(gca,'LineWidth',lw,'FontSize',axisFS);
xlabel('trial','FontSize',labelFS);
ylabel('confidence in others','FontSize',labelFS);
set(gca,'XTick',linspace(2,20,7));
set(gca,'YTick',[.5:.1:1]);
ylim([.5 1]);
xlim([0 21]);
box off;
print('-djpeg','-r300',['Figures/Figure7A_confidence']);