%% NetworkVisualizationAnalysis
% This script creates visualizations comparing real and fitted network characteristics,
% including bias analysis, border effects, and distance-based weight distributions.
%
% Input required:
% - BiasReal.mat: Contains bias measures from real network
% - BiasFitted.mat: Contains bias measures from fitted network
% - ConfiniFitted.mat: Contains analyzed fitted patterns by distance and borders
% - ConfiniReal.mat: Contains analyzed real patterns by distance and borders
%
% Outputs:
% - Figure 1: Bias comparison and border region analysis
% - Figure 2: Distance-based weight distribution analysis

%% Initialize Environment
clc;clearvars;close all

%% Load Analysis Results
% Load bias and network analysis results
load BiasReal.mat        % Contains bias measures for real network
load BiasFitted.mat      % Contains bias measures for fitted network
load ConfiniFitted.mat   % Contains fitted network analysis results
load ConfiniReal.mat     % Contains real network analysis results

% Find common countries between real and fitted analyses
[paesi,ia,ib]=intersect(paesi,paesiFitted);

%% Create Bias and Border Analysis Figure
figure

% Plot 1: Real Network Bias
subplot(2,2,1)
bar(bias)
title('Real Network')
ylabel('Within/Total % of Weights')
legend({'Country','Sector','Country-Sector'},'Location','best')
axis tight
grid on
set(findall(gcf,'-property','FontSize'),'FontSize',12)
xticks(1:length(Tempi))
xticklabels(Tempi)
colororder('sail')

% Plot 2: Fitted Network Bias
subplot(2,2,2)
bar(biasFitted)
title('Fitted Network')
ylabel('Within/Total % of Weights')
legend({'Country','Sector','Country-Sector'},'Location','best')
axis tight
grid on
set(findall(gcf,'-property','FontSize'),'FontSize',12)
xticks(1:length(Tempi))
xticklabels(Tempi)
colororder('sail')

% Plot 3: Real Network Border Analysis
subplot(2,2,3)
bar([log(varBP(ia)),log(varBP2(ia))])
axis tight
grid on
ylabel('Log Average Weights')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
xticks(1:size(varBP,1))
xticklabels(paesi)
title('Log Average Weights of Border NUTS-3 Inside vs Outside Country Real Network')
legend({'Inside','Outside'})
colororder('sail')

% Plot 4: Fitted Network Border Analysis
subplot(2,2,4)
bar([log(varBPfitted(ib)),log(varBP2fitted(ib))])
axis tight
grid on
ylabel('Log Average Weights')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
xticks(1:size(varBP,1))
xticklabels(paesi)
title('Log Average Weights of Border NUTS-3 Inside vs Outside Country Fitted Network')
legend({'Inside','Outside'})
colororder('sail')

% Statistical Tests for Border Effects
% Kolmogorov-Smirnov tests for distribution differences
[h p] = kstest2(log(varBP(ia)),log(varBP2(ia)))        % Real network
[h,p] = kstest2(log(varBPfitted(ib)),log(varBP2fitted(ib)))  % Fitted network

%% Create Distance-Based Weight Distribution Figure
% Initialize figure with specific properties
figure('Position', [100, 100, 1200, 800], 'Color', [0.95, 0.95, 0.95]);

% Define color scheme for visualization
colors = [0.2, 0.7, 0.8;   % Light blue
          0.8, 0.4, 0.2;   % Orange
          0.2, 0.7, 0.8;   % Light blue
          0.8, 0.4, 0.2];  % Orange

% Define common text properties
commonProps = {'FontWeight', 'bold', 'FontSize', 12};

% Create subplots for each distance quartile
for i = 1:4
    subplot(2, 2, i)
    
    % Compute data for current distance quartile
    data = eval(['[sum(Amount', num2str(i), 'b),sum(Amount', num2str(i), ...
        'a),sum(Amount', num2str(i), 'bfitted),sum(Amount', num2str(i), 'afitted)]']);
    
    % Create and customize bar plot
    b = bar(data, 'FaceColor', 'flat');
    b.CData = colors;
    
    % Configure x-axis
    xticks(1:4)
    xticklabels({'Same Country Real', 'Diff. Country Real', ...
        'Same Country Fitted', 'Diff. Country Fitted'})
    set(gca, 'XTickLabelRotation', 45)
    
    % Set y-axis limits with headroom
    ylim([0, max(data) * 1.2])
    
    % Add value labels on bars
    text(1:4, data, num2str(data', '%.0f'), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 10, 'FontWeight', 'bold');
    
    % Label y-axis
    ylabel('Weights Sum', commonProps{:})
    
    % Add subplot titles
    titles = {
        'Distance < 25-prctile', 
        '25-prctile < Distance < 50-prctile', 
        '50-prctile < Distance < 75-prctile', 
        'Distance > 75-prctile'
    };
    title(titles{i}, commonProps{:}, 'FontSize', 11);
    
    % Customize grid and appearance
    grid on
    set(gca, 'GridAlpha', 0.3, 'GridLineStyle', ':')
    box on
    set(findall(gca, '-property', 'FontSize'), 'FontSize', 11)
end

% Add overall title
sgtitle('Weight Distribution Across Distance Percentiles', ...
    'FontSize', 14, 'FontWeight', 'bold')

% Adjust subplot layout
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);