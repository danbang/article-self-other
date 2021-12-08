% Bang et al (2021) Neurocomputational mechanisms of confidence in self and
% other
%
% Reproduces Figure S4
%
% The simulated datasets (200 simulations x 21 subjects x 3 true models) and 
% the associated fits (200 simulations x 3 true models x 3 fitted models) 
% are not included in the repository due to overall file size (but will be
% shared upon request)
%
% The script instead loads pre-computed model metrics (WAIC and LOO) for
% each simulation x subject x true model x fitted model
%
% Dan Bang danbang.db@gmail.com 2021

%% -----------------------------------------------------------------------
%% PREPARATION

% fresh memory
clear; close all;

% Paths [change 'repoBase' according to local setup]
fs= filesep;
repoBase= [getDropbox(1),fs,'Ego',fs,'Matlab',fs,'ucl',fs,'self_other',fs,'Repository'];
dirDataBehaviour= [repoBase,fs,'Data',fs,'Task',fs,'Main'];
dirRecovery= [repoBase,fs,'Modelling',fs,'Recovery'];
dirFigures= [repoBase,fs,'Figures'];

% Specifications
modelz= {'S_B2T2P1','Q_B2T2P1', 'T_B2T2P1'};
n_sim= 200;

%% -----------------------------------------------------------------------
%% ANALYSIS

% Load pre-computed metrics
load([dirRecovery,fs,'precomputed_metrics.mat']);

% Loop through TRUE models
for i_true= 1:size(modelz,2);
        
    % Loop through iterations
    for i_sim= 1:n_sim
        
        
        % Loop through FIT models
        for i_fit= 1:size(modelz,2);
            
            % Log metrics  
            tmp_WAIC(i_fit)= mean(WAIC(i_sim,:,i_true,i_fit));
            tmp_LOO(i_fit)= mean(LOO(i_sim,:,i_true,i_fit));
        
        end    
        
        % Identify winner
        [v i]= min(tmp_WAIC);
        tmp= zeros(1,size(modelz,2));
        tmp(i)= 1;
        winner_WAIC(i_sim,:,i_true)= tmp;
        [v i]= min(tmp_LOO);
        tmp= zeros(1,size(modelz,2));
        tmp(i)= 1;
        winner_LOO(i_sim,:,i_true)= tmp;

    end
    
end

% Translate into fitted x true model probabilities
for i_true= 1:3;
   my_data(:,i_true)= sum(winner_WAIC(:,:,i_true));
end

%% -----------------------------------------------------------------------
%% VISUALISATION

%% FIGURE: CONFUSION
% Data
map= my_data;
% normalise
for c= 1:size(modelz,2);
        map(:,c)= map(:,c)./sum(map(:,c));
end
map(isnan(map))= 0;
% Plot
figure('color',[1 1 1]);
colormap('parula');
imagesc(map);
% add percentages onto plot
for r= 1:3;
    for c= 1:3;
        my_string= sprintf('%.2f', round(map(r,c)*100)/100);
        text(c-.17,r,my_string,'FontSize',20);
    end
end
% Tidy up
set(gca,'YTick',1:3,'YTickLabel',{'S','Q','T'});
set(gca,'XTick',1:3,'XTickLabel',{'S','Q','T'});
xlabel('true model');
ylabel('fitted model');
set(gca,'FontSize',24,'LineWidth',4);
c=colorbar;
set(c,'YTick',0:.2:1);
set(c,'LineWidth',4);
caxis([0 1]);
axis square
title('confusion: P(fit | true)','FontWeight','normal');
print('-djpeg','-r300',[dirFigures,fs,'Figure-S5-confusion']);


%% FIGURE: INVERSION
% Data
map= my_data;
for r= 1:size(modelz,2)
    map(r,:)= map(r,:)./sum(map(r,:));
end
map(isnan(map))= 0;
% Plot
figure('color',[1 1 1]);
colormap('parula');
imagesc(map);
% add percentages onto plot
for r= 1:3;
    for c= 1:3;
        my_string= sprintf('%.2f', round(map(r,c)*100)/100);
        text(c-.17,r,my_string,'FontSize',20);
    end
end
% Tidy up
set(gca,'YTick',1:3,'YTickLabel',{'S','Q','T'});
set(gca,'XTick',1:3,'XTickLabel',{'S','Q','T'});
xlabel('true model');
ylabel('fitted model');
set(gca,'FontSize',24,'LineWidth',4);
c=colorbar;
set(c,'YTick',0:.2:1);
set(c,'LineWidth',4);
caxis([0 1]);
axis square;
title('inversion: P(true | fit)','FontWeight','normal');
print('-djpeg','-r300',[dirFigures,fs,'Figure-S5-confusion']);