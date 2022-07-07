% Bang et al (2021) Neurocomputational mechanisms of confidence in self and
% other
%
% Reproduces Figure S2
%
% Dan Bang danbang.db@gmail.com 2021

%% -----------------------------------------------------------------------
%% PREPARATION

% fresh memory
clear; close all;

% Paths [change 'repoBase' according to local setup]
fs= filesep;
repoBase= [getDropbox(1),fs,'Ego',fs,'Matlab',fs,'ucl',fs,'self_other',fs,'Repository'];
dirDataBehaviour= [repoBase,fs,'Data',fs,'Behaviour',fs,'Scan',fs,'Task'];
dirModelFits= [repoBase,fs,'Modelling',fs,'Stan_fits'];
dirFigures= [repoBase,fs,'Figures'];

% Models
modelz= {'S_B1T1P0','S_B1T2P0','S_B2T1P0','S_B2T2P0', ...
         'S_B1T1P1','S_B1T2P1','S_B2T1P1','S_B2T2P1_U', ...
         'Q_B1T1P0','Q_B1T2P0','Q_B2T1P0','Q_B2T2P0', ...
         'Q_B1T1P1','Q_B1T2P1','Q_B2T1P1','Q_B2T2P1_U', ...
         'T_B1T1P0','T_B1T2P0','T_B2T1P0','T_B2T2P0', ...
         'T_B1T1P1','T_B1T2P1','T_B2T1P1','T_B2T2P1_U'};
modelNamez= {'S-B1T1P0','S-B1T2P0','S-B2T1P0','S-B2T2P0', ...
         'S-B1T1P1','S-B1T2P1','S-B2T1P1','S-B2T2P1', ...
         'Q-B1T1P0','Q-B1T2P0','Q-B2T1P0','Q-B2T2P0', ...
         'Q-B1T1P1','Q-B1T2P1','Q-B2T1P1','Q-B2T2P1', ...
         'T-B1T1P0','T-B1T2P0','T-B2T1P0','T-B2T2P0', ...
         'T-B1T1P1','T-B1T2P1','T-B2T1P1','T-B2T2P1'};
method= 'vb';
Nreruns= 4;

%% -----------------------------------------------------------------------
%% ANALYSIS

% Loop through models
for i_m=1:size(modelz,2);

        % Initialise variables        
        tmp_metric_WAIC= [];
        tmp_metric_LOO= [];
        % Loop through reruns
        for i_r= 1:Nreruns
            % load model fit
            load([dirModelFits,fs,'model_',modelz{i_m},'_seed_',num2str(i_r),'_vb.mat']);
            % Log metrics 
            tmp_metric_WAIC= [tmp_metric_WAIC; subjectWAIC'];
            tmp_metric_LOO= [tmp_metric_LOO; subjectLOO'];
        end
        % Average across reruns
        metric_WAIC(i_m,:)= mean(tmp_metric_WAIC,1);
        metric_LOO(i_m,:)= mean(tmp_metric_LOO,1);
        
end

%% -----------------------------------------------------------------------
%% VISUALISATION

%% WAIC BAR PLOT
% data
my_data= metric_WAIC;
n= size(my_data,2);
m= size(my_data,1);
mu= mean(my_data');
sem= std(my_data')./sqrt(n);
% Plot
figure('color',[1 1 1]);
% Bars
barh(mu,'FaceColor',[.6 1 .6],'LineWidth',2); hold on;
% Overlay subjects
for i=1:m; plot([mu(i)-sem(i) mu(i)+sem(i)],[i i],'k-','LineWidth',2); hold on; end;
for i=1:m; randz= (unifrnd(0,1,1,n)-.5)/2; plot(my_data(i,:),i+randz,'k.','MarkerSize',10); hold on; end;
% Tidy up
set(gca,'ydir','reverse');
set(gca,'yaxislocation','left');
set(gca,'YTick',1:m,'YTickLabel',modelNamez);
xlabel('WAIC');
set(gca,'FontSize',12,'LineWidth',2);
ylim([0 m+1]);
box off;
print('-djpeg','-r300',['Figures',fs,'Figure-S2A-waic']);

%% WAIC HEATMAP
% Data
my_data= metric_WAIC;
% Create matrix for heatmap
j=0;
for c= 1:3;
    for r= 1:8
        j= j+1;
        map(r,c)= mean(my_data(j,:));
    end
end
map = rsa.util.rankTransform(map); % Rank transform (requires RSA toolbox)
% Plot
figure('color',[1 1 1]);
colormap('parula');
imagesc(map);
% Tidy up
set(gca,'YTick',1:8,'YTickLabel',{'B1T1P0','B1T2P0','B2T1P0','B2T2P0','B1T1P1','B1T2P1','B2T1P1','B2T2P1'});
set(gca,'XTick',[1:3],'XTickLabel',{'S','Q','T'});
set(gca,'FontSize',24,'LineWidth',4);
c=colorbar;
set(c,'LineWidth',4);
set(c,'YTick',[1 24],'YTickLabel',{'min','max'});
title('WAIC [rank-transformed]','FontWeight','normal');
print('-djpeg','-r300',['Figures',fs,'Figure-S2B-waic']);

%% WAIC BAR PLOT
% data
my_data= metric_LOO;
n= size(my_data,2);
m= size(my_data,1);
mu= mean(my_data');
sem= std(my_data')./sqrt(n);
% Plot
figure('color',[1 1 1]);
% Bars
barh(mu,'FaceColor',[.6 1 .6],'LineWidth',2); hold on;
% Overlay subjects
for i=1:m; plot([mu(i)-sem(i) mu(i)+sem(i)],[i i],'k-','LineWidth',2); hold on; end;
for i=1:m; randz= (unifrnd(0,1,1,n)-.5)/2; plot(my_data(i,:),i+randz,'k.','MarkerSize',10); hold on; end;
% Tidy up
set(gca,'ydir','reverse');
set(gca,'yaxislocation','left');
set(gca,'YTick',1:m,'YTickLabel',modelNamez);
xlabel('PSIS-LOO');
set(gca,'FontSize',12,'LineWidth',2);
ylim([0 m+1]);
box off;
print('-djpeg','-r300',['Figures',fs,'Figure-S2A-loo']);

%% WAIC HEATMAP
% Data
my_data= metric_LOO;
% Create matrix for heatmap
j=0;
for c= 1:3;
    for r= 1:8
        j= j+1;
        map(r,c)= mean(my_data(j,:));
    end
end
map = rsa.util.rankTransform(map); % Rank transform (requires RSA toolbox)
% Plot
figure('color',[1 1 1]);
colormap('parula');
imagesc(map);
% Tidy up
set(gca,'YTick',1:8,'YTickLabel',{'B1T1P0','B1T2P0','B2T1P0','B2T2P0','B1T1P1','B1T2P1','B2T1P1','B2T2P1'});
set(gca,'XTick',[1:3],'XTickLabel',{'S','Q','T'});
set(gca,'FontSize',24,'LineWidth',4);
c=colorbar;
set(c,'LineWidth',4);
set(c,'YTick',[1 24],'YTickLabel',{'min','max'});
title('PSIS-LOO [rank-transformed]','FontWeight','normal');
print('-djpeg','-r300',['Figures',fs,'Figure-S2B-loo']);

% write .csv source data file
csvheaders= {'loo_S_B1T1P0','loo_S_B1T2P0','loo_S_B2T1P0','loo_S_B2T2P0', ...
             'loo_S_B1T1P1','loo_S_B1T2P1','loo_S_B2T1P1','loo_S_B2T2P1_U', ...
             'loo_Q_B1T1P0','loo_Q_B1T2P0','loo_Q_B2T1P0','loo_Q_B2T2P0', ...
             'loo_Q_B1T1P1','loo_Q_B1T2P1','loo_Q_B2T1P1','loo_Q_B2T2P1_U', ...
             'loo_T_B1T1P0','loo_T_B1T2P0','loo_T_B2T1P0','loo_T_B2T2P0', ...
             'loo_T_B1T1P1','loo_T_B1T2P1','loo_T_B2T1P1','loo_T_B2T2P1_U', ...
             'waic_S_B1T1P0','waic_S_B1T2P0','waic_S_B2T1P0','waic_S_B2T2P0', ...
             'waic_S_B1T1P1','waic_S_B1T2P1','waic_S_B2T1P1','waic_S_B2T2P1_U', ...
             'waic_Q_B1T1P0','waic_Q_B1T2P0','waic_Q_B2T1P0','waic_Q_B2T2P0', ...
             'waic_Q_B1T1P1','waic_Q_B1T2P1','waic_Q_B2T1P1','waic_Q_B2T2P1_U', ...
             'waic_T_B1T1P0','waic_T_B1T2P0','waic_T_B2T1P0','waic_T_B2T2P0', ...
             'waic_T_B1T1P1','waic_T_B1T2P1','waic_T_B2T1P1','waic_T_B2T2P1_U'};
csvdata= [metric_LOO' metric_WAIC'];
csvwrite_with_headers('SourceData/FigureS2',csvdata,csvheaders);