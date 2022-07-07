% Bang et al (2021) Neurocomputational mechanisms of confidence in self and
% other
%
% Reproduces Figure S1
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
dirModelSumstats= [repoBase,fs,'Modelling',fs,'Models_sumstats'];
dirFigures= [repoBase,fs,'Figures'];

% Subjects
sbj_v= [22:27 29:32 34:42 44:45]; % subject numbers
k= 0; for j= sbj_v; k= k+1; sbj_name{k}= ['mriData_sub_',num2str(j)]; end % subject names

% Specify winning model for overlaying onto behavioural data
bestModel= 'T_B2T2P1';

%% -----------------------------------------------------------------------
%% ANALAYSIS

% initialise subject index
k = 0;

% loop through subjects
for i_sbj = 1:length(sbj_v)
        
    % update subject index
    k=k+1;
    
    %% Load and concatenate behavioural data
    % fresh memory
    clear data;
    % loop through blocks
    for i_blk = 1:3       
        % file name
        fname = [dirDataBehaviour,fs,sbj_name{k},'.mat'];
        % load file
        load(fname);
        % include
        incl_dec = ~isnan(DATA(i_blk).typeI.response);
        incl_gam = DATA(i_blk).typeII.side>0;
        include  = (incl_dec+incl_gam)==2;
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
    
    %% Analyse behavioural data (S: self; O: other)
    
    %% PERFORMANCE: AGGREGATE
    
    % Accuracy
    for kk= 1:3; S_accuracy_skill(i_sbj,kk) = nanmean(data.acc(data.self==1&data.kskill==kk)); end
    for kk= 1:3; O_accuracy_skill(i_sbj,kk) = nanmean(data.acc(data.self==0&data.kskill==kk)); end
    
    % RT
    for kk= 1:3; S_reactiontime_skill(i_sbj,kk) = nanmean(data.rt1(data.self==1&data.kskill==kk)); end
    for kk= 1:3; O_reactiontime_skill(i_sbj,kk) = nanmean(data.rt1(data.self==0&data.kskill==kk)); end

    %% PERFORMANCE: GLM
    
    % Accuracy
    % Self-trials
    indx = find(data.self==1);
    rt1  = zscore(log(data.rt1(indx)));
    gamb = data.gamble(indx);
    cacc = data.acc(indx);
    cohz = zscore(data.theta(indx));
    valz = zscore(data.rg(indx)-data.rs(indx));
    skil = (data.kskill(indx)-2)/2;
    X = [cohz; ...
         skil; ...
        ]';
    Y=cacc';
    B=glmfit(X,Y,'binomial','link','logit');
    betaS_acc(i_sbj,:) = B(2:end);   
    
    % Reaction time
    % Self-trials
    indx = find(data.self==1);
    rt1  = zscore(log(data.rt1(indx)));
    gamb = data.gamble(indx);
    cacc = data.acc(indx);
    cohz = zscore(data.theta(indx));
    valz = zscore(data.rg(indx)-data.rs(indx));
    skil = (data.kskill(indx)-2)/2;
    X = [cohz; ...
         skil; ...
        ]';
    Y=rt1';
    B=glmfit(X,Y);
    betaS_rt1(i_sbj,:) = B(2:end);   
    
    %% GAMBLE: AGGREGATE
    
    % Gamble: others' ability
    for kk= 1:3; S_gamble_skill(i_sbj,kk) = nanmean(data.gamble(data.self==1&data.kskill==kk)); end
    for kk= 1:3; O_gamble_skill(i_sbj,kk) = nanmean(data.gamble(data.self==0&data.kskill==kk)); end

    % Gamble: overall
    gamMean(i_sbj,:) = [nanmean(data.gamble(data.self==1)) nanmean(data.gamble(data.self==0))];
        
    % Gamble: coherence
    theta_med = median(data.theta);
    theta_bin = data.theta>=theta_med;
    gamCohS(i_sbj,1) = nanmean(data.gamble(data.self==1&theta_bin==0));
    gamCohS(i_sbj,2) = nanmean(data.gamble(data.self==1&theta_bin==1));
    gamCohO(i_sbj,1) = nanmean(data.gamble(data.self==0&theta_bin==0));
    gamCohO(i_sbj,2) = nanmean(data.gamble(data.self==0&theta_bin==1));
    
    % Gamble: value
    val= data.rg-data.rs;
    val_med = median(val);
    val_bin = val>=val_med;
    gamValS(i_sbj,1) = nanmean(data.gamble(data.self==1&val_bin==0));
    gamValS(i_sbj,2) = nanmean(data.gamble(data.self==1&val_bin==1));
    gamValO(i_sbj,1) = nanmean(data.gamble(data.self==0&val_bin==0));
    gamValO(i_sbj,2) = nanmean(data.gamble(data.self==0&val_bin==1));
    
    % Gamble: accuracy
    indx = find(data.self==1);
    cacc = data.acc(indx);
    pacc = [0 cacc(1:end-1)];
    gamb = data.gamble(indx);
    gamAccS(i_sbj,1) = nanmean(gamb(pacc==0));
    gamAccS(i_sbj,2) = nanmean(gamb(pacc==1));
    indx = find(data.self==0);
    cacc = data.acc(indx);
    pacc = [0 cacc(1:end-1)];
    gamb = data.gamble(indx);
    gamAccO(i_sbj,1) = nanmean(gamb(pacc==0));
    gamAccO(i_sbj,2) = nanmean(gamb(pacc==1));
    
    %% GAMBLE: GLM
    
    % Both trial types
    indx = 1:length(data.self);
    gamb = data.gamble(indx);
    type = data.self(indx)-.5;
    cacc = data.acc(indx);
    cohz = zscore(data.theta(indx));
    valz = zscore(data.rg(indx)-data.rs(indx));
    skil = (data.kskill(indx)-2)/2;
    X = [type; ... 
         cohz; ...
         valz; ...
         skil; ...
         cohz.*type; ...
         valz.*type; ...
         skil.*type; ...
        ]';
    Y=gamb';
    B=glmfit(X,Y,'binomial','link','logit');
    betaSO(i_sbj,:) = B(2:end);  
    
    % Self-trials
    indx = find(data.self==1);
    gamb = data.gamble(indx);
    cacc = data.acc(indx);
    cohz = zscore(data.theta(indx));
    valz = zscore(data.rg(indx)-data.rs(indx));
    skil = (data.kskill(indx)-2)/2;
    X = [cohz; ...
         valz; ...
         skil; ...
        ]';
    Y=gamb';
    B=glmfit(X,Y,'binomial','link','logit');
    betaS(i_sbj,:) = B(2:end);   
    
    % Other-trials
    indx = find(data.self==0);
    gamb = data.gamble(indx);
    cacc = data.acc(indx);
    cohz = zscore(data.theta(indx));
    valz = zscore(data.rg(indx)-data.rs(indx));
    skil = (data.kskill(indx)-2)/2;
    X = [cohz; ...
         valz; ...
         skil; ...
        ]';
    Y=gamb';
    B=glmfit(X,Y,'binomial','link','logit');
    betaO(i_sbj,:) = B(2:end);
    
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
ms= 20;
axisFS= 24;
labelFS= 36;
jitter=.5;

% FIGURE: GLM FOR BOTH TRIAL TYPES 
% Create figure
figure('color',[1 1 1]);
hold on;
% Reference line
plot([0 8],[0 0],'k-','LineWidth',lw);
% Group mean and SEM
steplotmulticoloffalpha(betaSO,repmat([.5 .5 .5],7,1),lw,1,alphaz);
% Overlay subjects
for i=1:7;
    y= betaSO(:,i);
    x= violaPoints(i,betaSO(:,i),jitter);
    scatter(x,y,ms,dcol,'filled','MarkerFaceAlpha',.4);
end
% Tidy up
xlim([0 8]);
ylim([-3.5 3.5]); 
set(gca,'XTick',1:7,'XTickLabel',{'T','C','R','A','TxC','TxR','TxA'});
set(gca,'YTick',-3:3);
set(gca,'LineWidth',lw,'FontSize',axisFS);
ylabel('beta [arb. units]','FontSize',labelFS);
xlabel('predictor','FontSize',labelFS);
print('-djpeg','-r300',['Figures/FigureS1A']);

% FIGURE: GLM SEPARATELY FOR EACH TRIAL TYPE
% Create figure
figure('color',[1 1 1]);
hold on;
% Reference line
plot([0 7],[0 0],'k-','LineWidth',lw);
% Group mean and SEM
% Self
steplotmulticoloffalpha(betaS,repmat(scol,3,1),lw,1,alphaz);
% Other
steplotmulticoloffalpha(betaO,repmat(ocol,3,1),lw,4,alphaz);
% Overlay subjects
% Self
for i=1:3;
    y= betaS(:,i);
    x= violaPoints(i,betaS(:,i),jitter);
    scatter(x,y,ms,dcol,'filled','MarkerFaceAlpha',.4);
end
% Other
for i=1:3;
    y= betaO(:,i);
    x= violaPoints(i+3,betaO(:,i),jitter);
    scatter(x,y,ms,dcol,'filled','MarkerFaceAlpha',.4);
end
% Tidy up
xlim([0 7]);
ylim([-3.5 3.5]); 
set(gca,'XTick',1:6,'XTickLabel',{'C','R','A'});
set(gca,'YTick',-3:3);
set(gca,'LineWidth',lw,'FontSize',axisFS);
ylabel('beta [arb. units]','FontSize',labelFS);
xlabel('predictor','FontSize',labelFS);
print('-djpeg','-r300',['Figures/FigureS1B']);

% write .csv source data file
csvheaders= {'panel_a_T','panel_a_C','panel_a_R','panel_a_A','panel_a_TxC','panel_a_TxR','panel_a_TxA' ...
             'panel_b_self_C','panel_b_self_R','panel_b_self_A','panel_b_other_C','panel_b_other_R','panel_b_other_A'};
csvdata= [betaSO betaS betaO];
csvwrite_with_headers('SourceData/FigureS1',csvdata,csvheaders);