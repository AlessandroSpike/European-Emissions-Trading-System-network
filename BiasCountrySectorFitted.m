% Gravity Model Network Bias Analysis
% This script processes gravity model outputs and calculates network bias metrics
% for country, sector, and country-sector combinations

clc;clearvars;close all

%% Add required paths
addpath('Initialdata')  % Directory containing input data
%% Load and process gravity model output data
load('OutputGravity.mat','iso_d','iso_o','fitted','year')
% Convert country codes to string arrays, removing header rows
iso_d = string(iso_d.iso_d(2:end));    % Destination country codes
iso_o = string(iso_o.iso_o(2:end));    % Origin country codes
fitted = fitted.xfittedvalues(2:end);   % Fitted values from gravity model
anni = year.year(2:end);               % Years

% Load partition data for bias calculations
load('AdjMatrixTime.mat','countrypartition',...
    'countrysectorpartition','sectorpartition','Tempi','C')
Call = C;  % Store complete list of countries/entities

%% Initialize bias matrices
% bias: stores weighted bias calculations
% biasBin: stores binary (unweighted) bias calculations
bias = zeros(length(Tempi),3);     % [country, sector, country-sector]
biasBin = zeros(length(Tempi),3);  % Binary version

%% Main processing loop
for t = 1:length(Tempi)
    % Extract data for current time period
    tengo = anni==(Tempi(t));
    fitted_t = fitted(tengo);      % Fitted values for current year
    iso_ot = iso_o(tengo);        % Origin countries
    isodt = iso_d(tengo);         % Destination countries
    
    % Create adjacency matrix from gravity model outputs
    Nomi = [iso_ot,isodt];
    [C,~,ib] = unique(Nomi, 'stable');
    ib = reshape(ib, size(Nomi));
    Adj = full(sparse(ib(:,1),ib(:,2),fitted_t,length(C),length(C)));
    
    % Remove inactive nodes (no connections)
    tolgo = find(sum(Adj)+sum(Adj')==0);
    AdjN = Adj;
    AdjN(tolgo,:) = [];
    AdjN(:,tolgo) = [];
    C(tolgo) = [];
    
    % Process country-level bias
    [countrypartition_t,ia,ib] = intersect(Call,C);
    countrypartition2 = countrypartition(ia);
    Adjt = AdjN(ib,ib);
    AdjBin = Adjt;
    AdjBin(AdjBin>0) = 1;
    
    % Calculate country bias (weighted and binary)
    [G,~] = findgroups(countrypartition2);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1),Adjt,G);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1),AdjCountryInstallation',G);
    bias(t,1) = 100*sum(diag(AdjCountryInstallation))/sum(sum(AdjCountryInstallation));
    
    % Binary country bias
    [G,~] = findgroups(countrypartition2);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1),AdjBin,G);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1),AdjCountryInstallation',G);
    biasBin(t,1) = 100*sum(diag(AdjCountryInstallation))/sum(sum(AdjCountryInstallation));
    
    % Process sector-level bias
    [sectorpartition_t,ia,ib] = intersect(Call,C);
    sectorpartition2 = sectorpartition(ia);
    Adjt = AdjN(ib,ib);
    AdjBin = Adjt;
    AdjBin(AdjBin>0) = 1;
    
    % Calculate sector bias (weighted and binary)
    [G,~] = findgroups(sectorpartition2);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1),Adjt,G);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1),AdjCountryInstallation',G);
    bias(t,2) = 100*sum(diag(AdjCountryInstallation))/sum(sum(AdjCountryInstallation));
    
    % Binary sector bias
    [G,~] = findgroups(sectorpartition2);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1),AdjBin,G);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1),AdjCountryInstallation',G);
    biasBin(t,2) = 100*sum(diag(AdjCountryInstallation))/sum(sum(AdjCountryInstallation));
    
    % Process country-sector combined bias
    [countrysectorpartition_t,ia,ib] = intersect(Call,C);
    countrysectorpartition2 = countrysectorpartition(ia);
    Adjt = AdjN(ib,ib);
    AdjBin = Adjt;
    AdjBin(AdjBin>0) = 1;
    
    % Calculate country-sector combined bias (weighted and binary)
    [G,~] = findgroups(countrysectorpartition2);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1),Adjt,G);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1),AdjCountryInstallation',G);
    bias(t,3) = 100*sum(diag(AdjCountryInstallation))/sum(sum(AdjCountryInstallation));
    
    % Binary country-sector bias
    [G,ID] = findgroups(countrysectorpartition2);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1),AdjBin,G);
    AdjCountryInstallation = splitapply(@(x)sum(x, 1),AdjCountryInstallation',G);
    biasBin(t,3) = 100*sum(diag(AdjCountryInstallation))/sum(sum(AdjCountryInstallation));
end

%% Save results
biasFitted = bias;
save('BiasFitted.mat','Tempi','biasFitted');