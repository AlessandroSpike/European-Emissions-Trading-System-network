clc; clearvars; close all

% ============================
% Script Description
% ============================
% This script analyzes transportation and cultural distance data across
% different regions over time. It uses optimal transport methods to assess
% community interactions based on given data. The script outputs plots
% showing average cultural distances within and between communities.
%
% ============================
% Inputs:
% - 'OutputGravity.mat': Contains:
%     - iso_d: Destination identifiers.
%     - iso_o: Origin identifiers.
%     - amount: Amount of transport or interaction.
%     - fitted: Fitted values from a gravity model.
%     - year: Corresponding years for the data.
%
% - 'AdjMatrixTime.mat': Contains:
%     - Tempi: Time periods for the analysis.
%     - C: List of community identifiers.
%     - Distance: Matrix of distances (not used in this script).
%
% - 'Border_Distance_Info.mat': Contains:
%     - Nomi: Border information (not used in this script).
%     - DummyBorder: Dummy variables for border information (not used).
%
% - 'culturaldistance.mat': Contains:
%     - ctr: Country codes.
%     - Various columns of cultural distance data.
%
% ============================
% Outputs:
% - Figures showing the average cultural distances within and between
%   communities based on different factors.
% - Optional: Saves community aggregation data to 'CommunityAgg.mat'.

% Add necessary paths for functions and initial data
addpath("functions")
addpath("Initialdata")

%% Load Data
% Load output gravity model data from a .mat file
load('OutputGravity.mat', 'iso_d', 'iso_o', 'amount', 'fitted', 'year')
% Convert relevant variables to string and extract data excluding the first entry
iso_d = string(iso_d.iso_d(2:end)); % Destination identifiers as strings
iso_o = string(iso_o.iso_o(2:end)); % Origin identifiers as strings
amount = amount.amount(2:end); % Amount of transport or interaction
fitted = fitted.xfittedvalues(2:end); % Fitted values from the model
year = year.year(2:end); % Years corresponding to the data

% Load additional matrices required for calculations
load('AdjMatrixTime.mat', 'Tempi', 'C', 'Distance') % Temporal data
load('Border_Distance_Info', 'Nomi', 'DummyBorder'); % Border info (not used here)
NomiBorder = Nomi; % Store for potential future use
Call = C; % Community identifiers

%% Initialize Matrices for Temporal Data
% Prepare matrices to hold adjacency and fitted values over time
AdjT = zeros(length(Call), length(Call), length(Tempi)); % Adjacency matrices
FittedT = zeros(length(Call), length(Call), length(Tempi)); % Fitted values matrices

% Loop through each time period to construct temporal adjacency matrices
for t = 1:length(Tempi)
    % Identify indices corresponding to the current time period
    tengo = year == (Tempi(t));

    % Extract relevant data for the current time period
    amount_t = abs(amount(tengo)); % Transport amounts for the time period
    fitted_t = abs(fitted(tengo)); % Fitted values for the time period
    iso_ot = iso_o(tengo); % Origin identifiers for the time period
    isodt = iso_d(tengo); % Destination identifiers for the time period
    
    % Combine origin and destination into one array for unique identification
    Nomi = [iso_ot, isodt];
    [C, ~, ib] = unique(Nomi, 'stable'); % Unique community identifiers
    ib = reshape(ib, size(Nomi)); % Reshape for indexing
    
    % Create sparse adjacency and fitted matrices
    Adj = full(sparse(ib(:, 1), ib(:, 2), amount_t, length(C), length(C))); % Adjacency matrix
    Fitted = full(sparse(ib(:, 1), ib(:, 2), fitted_t, length(C), length(C))); % Fitted values matrix
    
    % Initialize large matrices for storing data
    n_large = length(Call); % Number of large nodes (communities)
    large_matrix = zeros(n_large); % Initialize adjacency matrix
    large_matrix2 = zeros(n_large); % Initialize fitted values matrix
    
    % Find indices of small nodes in the large nodes array
    [~, small_indices] = ismember(C, Call);
    
    % Create a logical index matrix for the small matrix elements in the large matrix
    [I, J] = meshgrid(small_indices, small_indices);
    
    % Assign values from small matrix to large matrix using logical indexing
    large_matrix(sub2ind(size(large_matrix), I(:), J(:))) = Adj(:); % Fill adjacency values
    large_matrix2(sub2ind(size(large_matrix2), I(:), J(:))) = Fitted(:); % Fill fitted values

    % Store the adjacency and fitted matrices for the current time period
    AdjT(:, :, t) = large_matrix; % Store adjacency for time t
    FittedT(:, :, t) = large_matrix2; % Store fitted for time t
end

% Sum adjacency and fitted matrices across time to create total matrices
Amount = sum(AdjT, 3); % Total amount across all time periods
Fitted = sum(FittedT, 3); % Total fitted values across all time periods

%% Entropy Optimal Transport (OT)
gamma = 0.001; % Regularization parameter for optimal transport
Amount = Amount / sum(sum(Amount)); % Normalize amount matrix
Amount(isnan(Amount)) = 0; % Replace NaNs with zeros to avoid issues

% Calculate costs based on fitted values for optimal transport
C = max(max(Fitted)) ./ Fitted; % Cost matrix for transport

% Calculate in and out strengths for communities
InStr = sum(Amount); 
InStr = InStr / sum(InStr); % Normalize in-strength
OutStr = sum(Amount, 2); 
OutStr = OutStr / sum(OutStr); % Normalize out-strength

% Perform Sinkhorn Optimal Transport
[T, a, b, Err, disto] = Sinkhorn_OT(C, gamma, OutStr, InStr', 10^-5, 100); % Execute OT
[CiOT, QOT] = OTmodularity_dir(Amount, .3, T); % Compute community modularity

%% Cultural Distance Analysis
load culturaldistance.mat % Load cultural distance data
varianzaculturale = zeros(6, 1); % Preallocate variance array for 6 factors
figure % Create a figure for bar plots

% Loop through each cultural distance factor for analysis
for f = 1:6
    stati = culturaldistance.ctr; % Country codes
    pdi = nanmean(table2array(culturaldistance(:, f + 2)), 2); % Mean values for cultural distance
    varianzaculturale(f) = nanvar(pdi); % Variance calculation for the factor
    nome = culturaldistance.Properties.VariableNames{f + 2}; % Get the name of the current factor
    CallStati=char(Call);
    CallStati=string(CallStati(:,1:2));
    % Initialize a cell array to store assigned values for provinces
    province_values = zeros(size(CallStati));

    % Loop through each province to assign cultural distance values
    for i = 1:numel(CallStati)
        country_code = CallStati(i); % Extract country code from province name
        country_idx = find(ismember(stati, country_code)); % Find index in stati
        
        % Assign value if country code is found
        if ~isempty(country_idx)
            province_values(i) = pdi(country_idx(1)); % Assign mean cultural distance
        else
            province_values(i) = nan; % Assign NaN if not found
        end
    end

    % Calculate cultural distance matrix
    [I, J] = meshgrid(province_values, province_values); % Create grid for distance calculations
    DistanzaCulturale(:,:,1) = I; % Store cultural distance matrix
    DistanzaCulturale(:,:,2) = J; 
    DistanzaCulturale1 = (DistanzaCulturale(:,:,1) - DistanzaCulturale(:,:,2)).^2; % Squared differences
    DistanzaCulturaleTot(:,:,f) = DistanzaCulturale1 / varianzaculturale(f); % Normalize by variance
    
    % Sort and visualize results
    [commSort, pos] = sort(CiOT); % Sort communities by modularity
    DistanzaCulturaleS = DistanzaCulturale1(pos, :); % Sort cultural distances
    DistanzaCulturaleS = DistanzaCulturaleS(:, pos); % Sort again for proper indexing

    % Aggregate cultural distance by community
    AdjCountry = splitapply(@(x) nanmean(x, 1), DistanzaCulturaleS, commSort); % Aggregate within communities
    AdjCountry = splitapply(@(x) nanmean(x, 1), AdjCountry', commSort); % Aggregate between communities
    
    % Plot results for the current cultural distance factor
    subplot(2, 3, f)
    bar([diag(AdjCountry), mean(AdjCountry - diag(diag(AdjCountry)), 2)]) % Plot within and between community distances
    axis tight
    grid on
    xlabel('Community ID') % Label for x-axis
    legend('Within-Community', 'Between-Community') % Legend for the plot
    title(['Cultural Distance: ', nome]) % Title indicating the current factor
end

% Finalize cultural distance analysis by averaging across selected factors
DistanzaCulturaleTot1 = mean(DistanzaCulturaleTot(:,:,1:4), 3); % Average of first 4 factors
DistanzaCulturaleTot2 = mean(DistanzaCulturaleTot, 3); % Average of all factors
[commSort, pos] = sort(CiOT); % Sort communities again
DistanzaCulturaleTot1 = DistanzaCulturaleTot1(pos, :); % Sort average distances
DistanzaCulturaleTot1 = DistanzaCulturaleTot1(:, pos); % Sort again for proper indexing
AdjCountry = splitapply(@(x) nanmean(x, 1), DistanzaCulturaleTot1, commSort); % Aggregate by community
AdjCountry = splitapply(@(x) nanmean(x, 1), AdjCountry', commSort); % Aggregate across rows
DistanzaCulturaleTot2 = DistanzaCulturaleTot2(pos, :); % Sort total distances
DistanzaCulturaleTot2 = DistanzaCulturaleTot2(:, pos); % Sort again for proper indexing
AdjCountry2 = splitapply(@(x) nanmean(x, 1), DistanzaCulturaleTot2, commSort); % Final aggregation
AdjCountry2 = splitapply(@(x) nanmean(x, 1), AdjCountry2', commSort); % Final aggregation across rows

% Create final figures to display average cultural distances
figure
subplot(1, 2, 1)
bar([diag(AdjCountry), mean(AdjCountry - diag(diag(AdjCountry)), 2)]) % Bar plot for 4 factors
axis tight
grid on
xlabel('Community ID')
legend('Within-Community', 'Between-Community')
title(['Average Cultural Distance 4 Factors'])

subplot(1, 2, 2)
bar([diag(AdjCountry2), mean(AdjCountry2 - diag(diag(AdjCountry2)), 2)]) % Bar plot for 6 factors
axis tight
grid on
xlabel('Community ID')
legend('Within-Community', 'Between-Community')
title(['Average Cultural Distance 6 Factors'])
