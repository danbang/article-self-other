% Bang et al (2021) Neurocomputational mechanisms of confidence in self and
% other
%
% Reproduces Figure S4
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

% Subjects
sbj_v= [22:27 29:32 34:42 44:45]; % subject numbers
k= 0; for j= sbj_v; k= k+1; sbj_name{k}= ['mriData_sub_',num2str(j)]; end % subject names

% Specify winning model for overlaying onto behavioural data
bestModel= 'T_B2T2P1';
Nreruns= 4;

%% -----------------------------------------------------------------------
%% ANALAYSIS

% initialise subject index
k = 0;

% initialise matrices
evolution.condition= NaN(length(sbj_v),120);
evolution.type= NaN(length(sbj_v),120);
evolution.empGamble= NaN(length(sbj_v),120);
evolution.simGamble= NaN(length(sbj_v),120);
evolution.simCon= NaN(length(sbj_v),120);
evolution.simNoise= NaN(length(sbj_v),120);

%% Extract data for group-level analysis
% loop through subjects
for i_sbj = 1:length(sbj_v);
        
    % update subject index
    k=k+1;
    
    %% Load and concatenate behavioural data
    % fresh memory
    clear data tmp;
    % loop through blocks
    for i_blk = 1:3       
        % file name
        fname = [dirDataBehaviour,fs,sbj_name{i_sbj},'.mat'];
        % load file
        load(fname);
        % include - modified to include all trials with additional vector
        % indicating which trials to keep
        incl_dec = ~isnan(DATA(i_blk).typeI.response);
        incl_gam = DATA(i_blk).typeII.side>0;
        incl_RT1 = ~isnan(DATA(i_blk).typeI.response_time);
        incl_RT2 = ~isnan(DATA(i_blk).typeII.response_time);
        include = ones(1,40)==1;
        tmp.include  = (incl_dec+incl_gam+incl_RT1+incl_RT2)==4;
        % load variables
        tmp.kskill = ones(1,sum(include)).*DATA(i_blk).trials.otherK;
        tmp.self = DATA(i_blk).trials.self(include);
        tmp.theta = DATA(i_blk).trials.theta(include);
        tmp.stmright = DATA(i_blk).trials.d(include)==1;
        tmp.decision = DATA(i_blk).typeI.response(include)-1;
        tmp.rs = DATA(i_blk).trials.rs(include);
        tmp.rg = DATA(i_blk).trials.rg(include);
        tmp.acc  = DATA(i_blk).typeI.correct(include);
        tmp.accO = DATA(i_blk).trials.otherCorrect(include);
        tmp.gamble = DATA(i_blk).typeII.response(include); 
        tmp.outcome = DATA(i_blk).typeII.r(include); 
        tmp.rt1  = DATA(i_blk).typeI.response_time(include); 
        tmp.rt1O = DATA(i_blk).trials.other_rt(include); 
        tmp.rt2  = DATA(i_blk).typeII.response_time(include); 
        % get data field names
        fn = fieldnames(tmp);
        % if first block, then initialise temporary storage structure
        if i_blk == 1; for i_field = 1:length(fn); eval(['data.',fn{i_field},'=[];']); end; end
        % add data to temporary storage structure
        for i_field = 1:length(fn); eval(['data.',fn{i_field},'=[data.',fn{i_field},' tmp.',fn{i_field},'];']); end
    end
    
    %% Load model data
    clear sim tmp;
    % loop through reruns
    for i_r= 1:Nreruns;
        load([dirModelFits,fs,'model_',bestModel,'_seed_',num2str(i_r),'_vb.mat']);
        tmp= subject_predict;
        fnames=fieldnames(tmp); 
        for i_f= 1:length(fnames); evalc(['sim.',fnames{i_f},'(i_r,:)=tmp.',fnames{i_f},'(i_sbj,:);']); end
    end
    % average over seeds
    fnames=fieldnames(tmp); 
    for i_f= 1:length(fnames); evalc(['sim.',fnames{i_f},'=mean(sim.',fnames{i_f},');']); end
    % re-arrange to match NaN trials between empirical and simulated data
    tmp_gamble= NaN(1,120);
    tmp_confidence= NaN(1,120);
    tmp_noise= NaN(1,120);
    data.trials= 1:120;
    for i_block= 1:3;
        c_indx= data.trials( ((i_block-1)*40+1): (i_block*40) );
        c_incl= data.include( ((i_block-1)*40+1): (i_block*40) );
        tmp_gamble(c_indx(c_incl==1))= sim.ppPgamble( ((i_block-1)*40+1): ((i_block*40)-(40-sum(c_incl))));
        tmp_confidence(c_indx(c_incl==1))= sim.ppCo( ((i_block-1)*40+1): ((i_block*40)-(40-sum(c_incl))));
        tmp_noise(c_indx(c_incl==1))= sim.ppk( ((i_block-1)*40+1): ((i_block*40)-(40-sum(c_incl))));
    end
    sim.ppPgamble= tmp_gamble;
    sim.ppCo= tmp_confidence;
    sim.ppk= tmp_noise;
    
    %% Save relevant data
    evolution.condition(i_sbj,:)= data.kskill;
    evolution.type(i_sbj,:)= data.self;
    evolution.empGamble(i_sbj,data.include==1)= data.gamble(data.include==1);
    evolution.simGamble(i_sbj,data.include==1)= sim.ppPgamble(data.include==1);
    evolution.simCon(i_sbj,data.include==1)= sim.ppCo(data.include==1);
    evolution.simNoise(i_sbj,data.include==1)= sim.ppk(data.include==1);
    
end

%% Group-level analysis: sliding window
% Evolution of p(gamble) by others' ability
windowSize= 3;
% loop through subjects
for i_sbj= 1:length(sbj_v)
    % loop through skill
    for i_skill= 1:3;
        my_empGamble= evolution.empGamble(i_sbj,evolution.condition(i_sbj,:)==i_skill & evolution.type(i_sbj,:)==0);
        my_simGamble= evolution.simGamble(i_sbj,evolution.condition(i_sbj,:)==i_skill & evolution.type(i_sbj,:)==0); 
        my_simCon= evolution.simCon(i_sbj,evolution.condition(i_sbj,:)==i_skill & evolution.type(i_sbj,:)==0); 
        my_simNoise= evolution.simNoise(i_sbj,evolution.condition(i_sbj,:)==i_skill & evolution.type(i_sbj,:)==0); 
        % moving average
        for t= windowSize:size(my_empGamble,2)
            S_my_empGamble(:,t-(windowSize-1))= nanmean(my_empGamble(:,(t+1-windowSize):t),2);
            S_my_simGamble(:,t-(windowSize-1))= nanmean(my_simGamble(:,(t+1-windowSize):t),2);
            S_my_simCon(:,t-(windowSize-1))= nanmean(my_simCon(:,(t+1-windowSize):t),2);
            S_my_simNoise(:,t-(windowSize-1))= nanmean(my_simNoise(:,(t+1-windowSize):t),2);
        end
        sliding.empGamble(i_sbj,:,i_skill)= S_my_empGamble;
        sliding.simGamble(i_sbj,:,i_skill)= S_my_simGamble;
        sliding.simCon(i_sbj,:,i_skill)= S_my_simCon;
        sliding.simNoise(i_sbj,:,i_skill)= S_my_simNoise;
    end
end

%% Group-level analysis: two halves
% Evolution of p(gamble) by others' ability
windowSize= 3;
% loop through subjects
for i_sbj= 1:length(sbj_v)
    % loop through skill
    for i_skill= 1:3;
        my_empGamble= evolution.empGamble(i_sbj,evolution.condition(i_sbj,:)==i_skill & evolution.type(i_sbj,:)==0);
        my_simGamble= evolution.simGamble(i_sbj,evolution.condition(i_sbj,:)==i_skill & evolution.type(i_sbj,:)==0); 
        halves.empGamble(i_sbj,1,i_skill)= nanmean(my_empGamble(1:10));
        halves.empGamble(i_sbj,2,i_skill)= nanmean(my_empGamble(11:20));
        halves.simGamble(i_sbj,1,i_skill)= nanmean(my_simGamble(1:10));
        halves.simGamble(i_sbj,2,i_skill)= nanmean(my_simGamble(11:20));
    end
end

%% -----------------------------------------------------------------------
%% VISUALISATION

%% General specifications
alphaz= .6;
lw= 4;
ms= 100;
axisFS= 20;
labelFS= 28;
jitter=.5;

%% P(GAMBLE) EVOLUTION
% Create figure
fig=figure('color',[1 1 1],'pos',[10 10 1200 300]);
fig.PaperPositionMode = 'auto';
hold on;
% Plot data
subplot(1,3,1);
i_skill= 1;
ci95fillplotmulticolalpha(squeeze(sliding.empGamble(:,:,i_skill)),lw,'-',[0 0 0]);
ci95fillplotmulticolalpha(squeeze(sliding.simGamble(:,:,i_skill)),lw,'-',[0 1 0]);
ylim([.3 .9]);
set(gca,'YTick',[.4 .6 .8],'XTick',[3 8 13 18],'XTickLabel',{'5','10','15','20'})
set(gca,'LineWidth',lw,'FontSize',axisFS);
ylabel('P(gamble)','FontSize',labelFS);
xlabel('trial','FontSize',labelFS);
title('low ability','FontWeight','normal','FontSize',labelFS);
subplot(1,3,2);
i_skill= 2;
ci95fillplotmulticolalpha(squeeze(sliding.empGamble(:,:,i_skill)),lw,'-',[0 0 0]);
ci95fillplotmulticolalpha(squeeze(sliding.simGamble(:,:,i_skill)),lw,'-',[0 1 0]);
ylim([.3 .9]);
set(gca,'YTick',[.4 .6 .8],'XTick',[3 8 13 18],'XTickLabel',{'5','10','15','20'})
set(gca,'LineWidth',lw,'FontSize',axisFS);
ylabel('P(gamble)','FontSize',labelFS);
xlabel('trial','FontSize',labelFS);
title('medium ability','FontWeight','normal','FontSize',labelFS);
subplot(1,3,3);
i_skill= 3;
ci95fillplotmulticolalpha(squeeze(sliding.empGamble(:,:,i_skill)),lw,'-',[0 0 0]);
ci95fillplotmulticolalpha(squeeze(sliding.simGamble(:,:,i_skill)),lw,'-',[0 1 0]);
ylim([.3 .9]);
set(gca,'YTick',[.4 .6 .8],'XTick',[3 8 13 18],'XTickLabel',{'5','10','15','20'})
set(gca,'LineWidth',lw,'FontSize',axisFS);
ylabel('P(gamble)','FontSize',labelFS);
xlabel('trial','FontSize',labelFS);
title('high ability','FontWeight','normal','FontSize',labelFS);
print('-djpeg','-r300',['Figures/Figure-S4A']);

%% NOISE ESTIMATE EVOLUTION
% Create figure
fig=figure('color',[1 1 1],'pos',[10 10 1200 300]);
fig.PaperPositionMode = 'auto';
hold on;
% Plot data
subplot(1,3,1);
i_skill= 1;
ci95fillplotmulticolalpha(squeeze(sliding.simNoise(:,:,i_skill)),lw,'-',[0 1 0]);
ylim([0 .4]);
set(gca,'YTick',0:.1:.4,'XTick',[3 8 13 18],'XTickLabel',{'5','10','15','20'})
set(gca,'LineWidth',lw,'FontSize',axisFS);
ylabel('other''s sensory noise','FontSize',labelFS);
xlabel('trial','FontSize',labelFS);
title('low ability','FontWeight','normal','FontSize',labelFS);
subplot(1,3,2);
i_skill= 2;
ci95fillplotmulticolalpha(squeeze(sliding.simNoise(:,:,i_skill)),lw,'-',[0 1 0]);
ylim([0 .4]);
set(gca,'YTick',0:.1:.4,'XTick',[3 8 13 18],'XTickLabel',{'5','10','15','20'})
set(gca,'LineWidth',lw,'FontSize',axisFS);
ylabel('other''s sensory noise','FontSize',labelFS);
xlabel('trial','FontSize',labelFS);
title('medium ability','FontWeight','normal','FontSize',labelFS);
subplot(1,3,3);
i_skill= 3;
ci95fillplotmulticolalpha(squeeze(sliding.simNoise(:,:,i_skill)),lw,'-',[0 1 0]);
ylim([0 .4]);
set(gca,'YTick',0:.1:.4,'XTick',[3 8 13 18],'XTickLabel',{'5','10','15','20'})
set(gca,'LineWidth',lw,'FontSize',axisFS);
ylabel('other''s sensory noise','FontSize',labelFS);
xlabel('trial','FontSize',labelFS);
title('high ability','FontWeight','normal','FontSize',labelFS);
print('-djpeg','-r300',['Figures/Figure-S4B']);

%% NOISE ESTIMATE EVOLUTION
% Create figure
fig=figure('color',[1 1 1],'pos',[10 10 1200 300]);
fig.PaperPositionMode = 'auto';
hold on;
% Plot data
subplot(1,3,1);
i_skill= 1;
ci95fillplotmulticolalpha(squeeze(sliding.simCon(:,:,i_skill)),lw,'-',[0 1 0]);
ylim([.6 .9]);
set(gca,'YTick',.6:.1:.9,'XTick',[3 8 13 18],'XTickLabel',{'5','10','15','20'})
set(gca,'LineWidth',lw,'FontSize',axisFS);
ylabel('confidence in other','FontSize',labelFS);
xlabel('trial','FontSize',labelFS);
title('low ability','FontWeight','normal','FontSize',labelFS);
subplot(1,3,2);
i_skill= 2;
ci95fillplotmulticolalpha(squeeze(sliding.simCon(:,:,i_skill)),lw,'-',[0 1 0]);
ylim([.6 .9]);
set(gca,'YTick',.6:.1:.9,'XTick',[3 8 13 18],'XTickLabel',{'5','10','15','20'})
set(gca,'LineWidth',lw,'FontSize',axisFS);
ylabel('confidence in other','FontSize',labelFS);
xlabel('trial','FontSize',labelFS);
title('medium ability','FontWeight','normal','FontSize',labelFS);
subplot(1,3,3);
i_skill= 3;
ci95fillplotmulticolalpha(squeeze(sliding.simCon(:,:,i_skill)),lw,'-',[0 1 0]);
ylim([.6 .9]);
set(gca,'YTick',.6:.1:.9,'XTick',[3 8 13 18],'XTickLabel',{'5','10','15','20'})
set(gca,'LineWidth',lw,'FontSize',axisFS);
ylabel('confidence in other','FontSize',labelFS);
xlabel('trial','FontSize',labelFS);
title('high ability','FontWeight','normal','FontSize',labelFS);
print('-djpeg','-r300',['Figures/Figure-S4C']);