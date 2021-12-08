% Bang et al (2021) Neurocomputational mechanisms of confidence in self and
% other
%
% Extracts summary statistics for each model for overlaying onto
% behavioural data
%
% Dan Bang danbang.db@gmail.com 2021

%% -----------------------------------------------------------------------
%% PREPARATION

% fresh memory
clc; clear;

% Paths [change 'repoBase' according to local setup]
fs= filesep;
repoBase= [getDropbox(1),fs,'Ego',fs,'Matlab',fs,'ucl',fs,'self_other',fs,'Repository'];
dirDataBehaviour= [repoBase,fs,'Data',fs,'Behaviour',fs,'Scan',fs,'Task'];
dirModelFits= [repoBase,fs,'Modelling',fs,'Stan_fits'];
dirModelOutput= [repoBase,fs,'Modelling',fs,'Models_sumstats'];
dirFigures= [repoBase,fs,'Figures'];

% Add paths
addpath('Utils');

% Models
modelz= {'S_B1T1P0','S_B1T2P0','S_B2T1P0','S_B2T2P0', ...
         'S_B1T1P1','S_B1T2P1','S_B2T1P1','S_B2T2P1', ...
         'Q_B1T1P0','Q_B1T2P0','Q_B2T1P0','Q_B2T2P0', ...
         'Q_B1T1P1','Q_B1T2P1','Q_B2T1P1','Q_B2T2P1', ...
         'T_B1T1P0','T_B1T2P0','T_B2T1P0','T_B2T2P0', ...
         'T_B1T1P1','T_B1T2P1','T_B2T1P1','T_B2T2P1'};
modelNamez= {'S-B1T1P0','S-B1T2P0','S-B2T1P0','S-B2T2P0', ...
         'S-B1T1P1','S-B1T2P1','S-B2T1P1','S-B2T2P1', ...
         'Q-B1T1P0','Q-B1T2P0','Q-B2T1P0','Q-B2T2P0', ...
         'Q-B1T1P1','Q-B1T2P1','Q-B2T1P1','Q-B2T2P1', ...
         'T-B1T1P0','T-B1T2P0','T-B2T1P0','T-B2T2P0', ...
         'T-B1T1P1','T-B1T2P1','T-B2T1P1','T-B2T2P1'};
method= 'vb';
Nreruns= 4;

% Subjects
sbj_v= [22:27 29:32 34:42 44:45]; % subject numbers
k= 0; for j= sbj_v; k= k+1; sbj_name{k}= ['mriData_sub_',num2str(j)]; end % subject names

%% -----------------------------------------------------------------------
%% ANALYSIS

% Loop through models
for i_m=1:size(modelz,2);
   
    k = k+1;
    
    % Loop through subjects
    for i_s= 1:size(sbj_v,2)
        
        % Load and parse empirical data
        load([dirDataBehaviour,filesep,sbj_name{i_s},'.mat']);
        concatenate;
        fnames=fieldnames(tmp); 
        for i_f= 1:length(fnames); evalc(['data.',fnames{i_f},'=tmp.',fnames{i_f},'(tmp.incl==1);']); end
        
        % Load and parse model data
        % Loop through reruns
        sim_pgamble = []; % p(gamble) samples for each trial = Nreruns x outsamp
        for i_r= 1:Nreruns
            load([dirModelFits,fs,'model_',modelz{i_m},'_seed_',num2str(i_r),'_vb.mat']);
            nonanz= subject_predict.lik(i_s,:)<1; % exclude NaN trials
            sim_pgamble= [sim_pgamble; squeeze(ppPgamble(:,i_s,nonanz))];
        end

        % loop through samples
        for kk= 1:size(sim_pgamble,1)
            % coherence
            theta_med = median(data.theta);
            theta_bin = data.theta>=theta_med;
            gamb = sim_pgamble(kk,:);
            gamCohS(kk,1) = nanmean(gamb(data.self==1&theta_bin==0));
            gamCohS(kk,2) = nanmean(gamb(data.self==1&theta_bin==1));
            gamCohO(kk,1) = nanmean(gamb(data.self==0&theta_bin==0));
            gamCohO(kk,2) = nanmean(gamb(data.self==0&theta_bin==1));
            % value
            val = data.rg-data.rs;
            val_med = median(val);
            val_bin = val>=val_med;
            gamb = sim_pgamble(kk,:);
            gamValS(kk,1) = nanmean(gamb(data.self==1&val_bin==0));
            gamValS(kk,2) = nanmean(gamb(data.self==1&val_bin==1));
            gamValO(kk,1) = nanmean(gamb(data.self==0&val_bin==0));
            gamValO(kk,2) = nanmean(gamb(data.self==0&val_bin==1));
            % skill
            gamb = sim_pgamble(kk,:);
            for jj= 1:3; gamSklS(kk,jj) = nanmean(gamb(data.self==1&data.kskill==jj)); end
            for jj= 1:3; gamSklO(kk,jj) = nanmean(gamb(data.self==0&data.kskill==jj)); end
            % accuracy
            indx = find(data.self==1);
            cacc = data.acc(indx);
            pacc = [0 cacc(1:end-1)];
            gamb = sim_pgamble(kk,indx);
            sim_sgam(kk,1) = nanmean(gamb(pacc==0));
            sim_sgam(kk,2) = nanmean(gamb(pacc==1));
            indx = find(data.self==0);
            cacc = data.acc(indx);
            pacc = [0 cacc(1:end-1)];
            gamb = sim_pgamble(kk,indx);
            sim_ogam(kk,1) = nanmean(gamb(pacc==0));
            sim_ogam(kk,2) = nanmean(gamb(pacc==1));
        end
        modz{i_m}.gamAccS(i_s,:)= mean(sim_sgam);
        modz{i_m}.gamAccO(i_s,:)= mean(sim_ogam);
        modz{i_m}.gamCohS(i_s,:)= mean(gamCohS);
        modz{i_m}.gamCohO(i_s,:)= mean(gamCohO);
        modz{i_m}.gamSklS(i_s,:)= mean(gamSklS);
        modz{i_m}.gamSklO(i_s,:)= mean(gamSklO);
        modz{i_m}.gamValS(i_s,:)= mean(gamValS);
        modz{i_m}.gamValO(i_s,:)= mean(gamValO);
        
    end
    
    my_model= modz{i_m};
    save(['Models_sumstats/',modelz{i_m},'.mat'],'my_model');

end