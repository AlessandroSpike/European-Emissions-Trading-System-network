% This MATLAB script analyzes and visualizes deviance data for different gravity models
% of carbon emission trading. It compares various model configurations and
% calculates confidence intervals for the difference in deviances.

% Input:
% - CSV files containing deviance data for different model configurations
% - Years range from 2005 to 2020

% Output:
% - Two figures: 
%   1. Deviance plots for different model configurations
%   2. Log of deviance difference with confidence intervals
% - Console output of difference in deviances and confidence intervals

clc;clearvars;close all

%% Add path to data files
addpath("Initialdata\GravityModelPerformance\")

%% Load data from CSV files
% Load and extract deviance data for various model configurations
DevinceBasic=readmatrix("db_deviance_only_standard_controls.csv");
DevinceBasic=DevinceBasic(2,2:end);

DevinceBasic_noEmiss=readmatrix("db_deviance_full_controls_no_emissions.csv");
DevinceBasic_noEmiss=DevinceBasic_noEmiss(2,2:end);

DevinceBasic_full_control=readmatrix("db_deviance_full_controls.csv");
DevinceBasic_full_control=DevinceBasic_full_control(2,2:end);

DeviancePaesi=readmatrix("db_deviance_full_controls_and_country.csv");
DeviancePaesi=DeviancePaesi(2,2:end);

DeviancePaesiSec=readmatrix("db_deviance.csv");
DeviancePaesiSec=DeviancePaesiSec(2,2:end);

DevianceKS=readmatrix("db_deviance_ks.csv");
DevianceKS=DevianceKS(2,2:end);

%% Calculate confidence bounds
% Given deviances and degrees of freedom
D1 = DeviancePaesiSec;  % Deviance for simpler model
D2 = DevianceKS;  % Deviance for more complex model

% Set confidence level (e.g., 95%)
alpha = 0.05;

% Calculate difference in deviance and degrees of freedom
delta_D = D1 - D2;
delta_df = ones(1,16);

% Calculate Chi-squared critical values for confidence bounds
chi2_low = chi2inv(1 - alpha/2, delta_df);  % Upper critical value
chi2_high = chi2inv(alpha/2, delta_df);     % Lower critical value

% Calculate confidence interval for the difference in deviances
lower_bound = delta_D ./ chi2_low;
upper_bound = delta_D ./ chi2_high;

% Display the results
disp(['Difference in Deviances: ', num2str(delta_D)]);
disp(['95% Confidence Interval for Difference in Deviances: [', num2str(lower_bound), ', ', num2str(upper_bound), ']']);

%% Initialize variables
Anni=2005:2020;

%% Create plots
% Plot deviance for different model configurations
figure
plot(Anni,DevinceBasic,'LineWidth',2)
hold on
plot(Anni,DevinceBasic_noEmiss,'LineWidth',2)
plot(Anni,DevinceBasic_full_control,'LineWidth',2)
plot(Anni,DeviancePaesi,'LineWidth',2)
plot(Anni,DeviancePaesiSec,'LineWidth',2)
xline(Anni(1),'--','Phase 1','LabelVerticalAlignment', 'middle')
xline(Anni(4),'--','Phase 2','LabelVerticalAlignment', 'middle')
xline(Anni(9),'--','Phase 3','LabelVerticalAlignment', 'middle')
axis tight
grid on
legend('Basic','Allocation','Emission','Country','Sector')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
title('Deviance')
colororder("meadow")

% Plot log of deviance difference with confidence intervals
figure
plot(Anni,log(delta_D),'LineWidth',2)
hold on
plot(Anni,log(lower_bound),'r--','LineWidth',1)
hold on
plot(Anni,log(upper_bound),'r--','LineWidth',1)
xline(Anni(1),'--','Phase 1','LabelVerticalAlignment', 'middle')
xline(Anni(4),'--','Phase 2','LabelVerticalAlignment', 'middle')
xline(Anni(9),'--','Phase 3','LabelVerticalAlignment', 'middle')
axis tight
grid on
set(findall(gcf,'-property','FontSize'),'FontSize',12)
title('Log of Deviance Difference')