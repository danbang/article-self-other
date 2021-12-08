% Bang et al (2021) Neurocomputational mechanisms of confidence in self and
% other
%
% Reproduces Figure S6A
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
modelz= {'S_B2T2P1_U','Q_B2T2P1_U','T_B2T2P1_U'};
modelNamez= {'S-B2T2P1','Q-B2T2P1','T-B2T2P1'};
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
mu= sum(my_data');
sem= std(my_data')./sqrt(n);
% Plot
figure('color',[1 1 1]);
% Bars
bar(mu,'FaceColor',[.6 1 .6],'LineWidth',2); hold on;
% Tidy up
set(gca,'FontSize',24,'LineWidth',4);
set(gca,'XTick',1:3,'XTickLabel',{'S','Q','ToM'});
ylabel('WAIC','FontSize',28);
xlabel('model + utility function','FontSize',28);
ylim([2400 2700]);
box off;
print('-djpeg','-r300',['Figures',fs,'Figure-S6A']);