clc;clearvars;close all

% Initialize workspace and clear previous data
%% Load required data
addpath('Initialdata')  % Add path to directory containing initial data files

%% Load adjacency matrix data
load AdjMatrixTime.mat  % Contains time-series adjacency matrices 'Adj' and timestamps 'Tempi'

%% Initialize bias matrices
% bias: stores weighted bias calculations
% biasBin: stores binary (unweighted) bias calculations
% Each matrix has 3 columns for different partition types:
% Col 1: Country bias
% Col 2: Sector bias
% Col 3: Country-Sector combined bias
bias = zeros(length(Tempi), 3);
biasBin = zeros(length(Tempi), 3);

%% Main processing loop - Calculate bias for each time point
for t = 1:length(Tempi)
    % Extract adjacency matrix for current time point
    Adjt = Adj(:,:,t);
    
    % Calculate country-level bias (weighted)
    [G,~] = findgroups(countrypartition);
    % Aggregate connections within country groups
    AdjCountryInstallation = splitapply(@(x)sum(x, 1), Adjt, G);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1), AdjCountryInstallation', G);
    % Calculate bias as percentage of within-group connections
    bias(t,1) = 100 * sum(diag(AdjCountryInstallation)) / sum(sum(AdjCountryInstallation));
    
    % Calculate sector-level bias (weighted)
    [G,~] = findgroups(sectorpartition);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1), Adjt, G);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1), AdjCountryInstallation', G);
    bias(t,2) = 100 * sum(diag(AdjCountryInstallation)) / sum(sum(AdjCountryInstallation));
    
    % Calculate country-sector combined bias (weighted)
    [G,ID] = findgroups(countrysectorpartition);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1), Adjt, G);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1), AdjCountryInstallation', G);
    bias(t,3) = 100 * sum(diag(AdjCountryInstallation)) / sum(sum(AdjCountryInstallation));
    
    % Create binary version of adjacency matrix (unweighted)
    AdjBin = Adjt;
    AdjBin(AdjBin > 0) = 1;
    
    % Repeat all bias calculations for binary matrix
    % Country-level bias (unweighted)
    [G,~] = findgroups(countrypartition);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1), AdjBin, G);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1), AdjCountryInstallation', G);
    biasBin(t,1) = 100 * sum(diag(AdjCountryInstallation)) / sum(sum(AdjCountryInstallation));
    
    % Sector-level bias (unweighted)
    [G,~] = findgroups(sectorpartition);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1), AdjBin, G);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1), AdjCountryInstallation', G);
    biasBin(t,2) = 100 * sum(diag(AdjCountryInstallation)) / sum(sum(AdjCountryInstallation));
    
    % Country-sector combined bias (unweighted)
    [G,ID] = findgroups(countrysectorpartition);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1), AdjBin, G);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1), AdjCountryInstallation', G);
    biasBin(t,3) = 100 * sum(diag(AdjCountryInstallation)) / sum(sum(AdjCountryInstallation));
end

%% Save results
save('BiasReal.mat', 'Tempi', 'bias');