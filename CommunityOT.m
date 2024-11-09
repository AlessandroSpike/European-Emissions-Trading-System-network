
% MATLAB Script for Network Analysis with Optimal Transport Methods
%
% Description:
% This script performs a network analysis over time using Optimal Transport
% (OT) methods and classical modularity approaches. It computes modularity
% and partition metrics for a set of temporal networks and compares them
% with community structures based on countries and sectors.
%
% Inputs:
% - OutputGravity.mat: Contains data on countries, fitted values, and amounts.
% - AdjMatrixTime.mat: Contains adjacency matrices and partition information.
% - CommunityAgg.mat: Provides aggregated community information.
%
% Outputs:
% - Plots and figures representing modularity and partition metrics over time.
% - Saved figures in 'PlotFig/' directory.

clc; 
clearvars; 
close all;

%% Add Necessary Paths
addpath("Initialdata"); % Path to input data files
addpath("functions");   % Path to custom functions used in the analysis

%% Load Data
% Load gravity model output data
load('OutputGravity.mat', 'iso_d', 'iso_o', 'fitted', 'amount', 'year');
iso_d = string(iso_d.iso_d(2:end)); % Source countries (exclude header)
iso_o = string(iso_o.iso_o(2:end)); % Destination countries (exclude header)
fitted = fitted.xfittedvalues(2:end); % Fitted values
amount = amount.amount(2:end); % Amounts
year = year.year(2:end); % Years

% Load adjacency matrix and partition data
load('AdjMatrixTime.mat', 'C', 'countrypartition', 'countrysectorpartition', 'sectorpartition');
Call = C; % Complete list of nodes

%% Initialize Variables
Tempi = 2005:2020; % Time period for analysis
gamma = [0.0001 0.001 0.01]; % Entropy regularization parameters for OT

% Containers for metrics
ModulOT = zeros(length(Tempi), length(gamma));
Modul = zeros(length(Tempi), 1);
Distanzona = zeros(length(Tempi), 1);
NummCompo = zeros(length(Tempi), 1);
CompoDistr = cell(length(Tempi), 1);
VariationInfoMetrics = zeros(length(Tempi), length(gamma));
MutualInfoMetrics = zeros(length(Tempi), length(gamma));
VariationInfoMetricsPaesi = zeros(length(Tempi), length(gamma) + 1);
MutualInfoMetricsPaesi = zeros(length(Tempi), length(gamma) + 1);
VariationInfoMetricsSettori = zeros(length(Tempi), length(gamma) + 1);
MutualInfoMetricsSettori = zeros(length(Tempi), length(gamma) + 1);
PartitionOT = cell(length(Tempi), length(gamma));
PartitioniOT1 = nan(length(Call), length(Tempi));
PartitioniOT2 = nan(length(Call), length(Tempi));
PartitioniOT3 = nan(length(Call), length(Tempi));
Partition = cell(length(Tempi), 1);
VariationInfoMetricsPaesiSettori = zeros(length(Tempi), length(gamma) + 1);
MutualInfoMetricsPaesiSettori = zeros(length(Tempi), length(gamma) + 1);

%% Main Loop Over Time
for t = 1:length(Tempi)
    disp(['Processing year: ', num2str(Tempi(t))]);
    tengo = year == (Tempi(t));

    % Extract data for the current year
    amount_t = abs(amount(tengo));
    fitted_t = abs(fitted(tengo));
    iso_ot = iso_o(tengo);
    isodt = iso_d(tengo);
    Nomi = [iso_ot, isodt];
    [C, ~, ib] = unique(Nomi, 'stable');
    ib = reshape(ib, size(Nomi));

    % Create adjacency matrices for actual and fitted amounts
    Adj = full(sparse(ib(:,1), ib(:,2), amount_t, length(C), length(C)));
    AdjFit = full(sparse(ib(:,1), ib(:,2), fitted_t, length(C), length(C)));

    % Initialize large matrices for the complete set of nodes
    n_large = length(Call);
    large_matrix = zeros(n_large); 
    large_matrixFit = zeros(n_large); 

    % Map small matrix indices to large matrix indices
    [Lia, small_indices] = ismember(C, Call); 
    indici = nan(size(Call));
    indici(small_indices) = small_indices;

    % Create a logical index matrix for mapping small matrix elements to the large matrix
    [I, J] = meshgrid(small_indices, small_indices);    
    % Assign values from the small adjacency matrix to the large adjacency matrix
    large_matrix(sub2ind(size(large_matrix), I(:), J(:))) = Adj(:);
    large_matrixFit(sub2ind(size(large_matrix), I(:), J(:))) = AdjFit(:); 
    AdjN = large_matrix; % Normalized adjacency matrix
    AdjNFit = large_matrixFit; % Fitted adjacency matrix

    % Remove nodes with zero in-degree and out-degree
    tolgo = find(sum(AdjN) + sum(AdjN') == 0);
    AdjN(tolgo, :) = [];
    AdjN(:, tolgo) = [];
    AdjNFit(tolgo, :) = [];
    AdjNFit(:, tolgo) = [];
    indici(tolgo) = [];
    C_nomi = C;

    % Update partitions after removing zero-degree nodes
    Country = countrypartition;
    Country(tolgo) = [];
    Sector = sectorpartition;
    Sector(tolgo) = [];
    CountrySector = countrysectorpartition;
    CountrySector(tolgo) = [];

    % Create directed graph and find connected components
    KK = digraph(AdjN);
    [bins, binsizes] = conncomp(KK);
    NummCompo(t, 1) = max(bins);
    CompoDistr{t, 1} = binsizes;
    [bins, binsizes] = conncomp(KK, 'Type', 'weak');
    NummCompo(t, 2) = max(bins);
    CompoDistr{t, 2} = binsizes;

    % Normalize adjacency matrix
    AdjN = AdjN / sum(sum(AdjN));
    AdjN(isnan(AdjN)) = 0;

    % Calculate cost matrix for OT
    C = 1 ./ AdjNFit;

    % Calculate in-strength and out-strength
    InStr = sum(AdjN); 
    InStr = InStr / sum(InStr);
    OutStr = sum(AdjN, 2); 
    OutStr = OutStr / sum(OutStr);

    % Perform Optimal Transport calculations for each gamma
    for k = 1:length(gamma)
        [T, a, b, Err, disto] = Sinkhorn_OT(C, gamma(k), OutStr, InStr', 10^-5, 100);
        Distanzona(t) = disto;

        % Calculate modularity using OT
        [CiOT, QOT] = OTmodularity_dir(AdjN, 1, T);
        ModulOT(t, k) = QOT;
        PartitionOT{t, k} = CiOT;

        % Store partition results for specific gamma values
        if k == 1
            PartitioniOT1(indici, t) = CiOT;
        elseif k == 2
            PartitioniOT2(indici, t) = CiOT;
        elseif k == 3
            PartitioniOT3(indici, t) = CiOT;
        end
    end
    
    % Calculate classical modularity
    [Ci, Q] = modularity_dir(AdjN, 1);
    Modul(t) = Q;
    Partition{t} = Ci;

    % Compare partitions and calculate variation and mutual information metrics
    for k = 1:length(gamma)
        [VIn, MIn] = partition_distance(Partition{t}, PartitionOT{t, k});
        VariationInfoMetrics(t, k) = VIn;
        MutualInfoMetrics(t, k) = MIn;

        [VIn, MIn] = partition_distance(Country, PartitionOT{t, k});
        VariationInfoMetricsPaesi(t, k) = VIn;
        MutualInfoMetricsPaesi(t, k) = MIn;

        [VIn, MIn] = partition_distance(Sector, PartitionOT{t, k});
        VariationInfoMetricsSettori(t, k) = VIn;
        MutualInfoMetricsSettori(t, k) = MIn;

        [VIn, MIn] = partition_distance(CountrySector, PartitionOT{t, k});
        VariationInfoMetricsPaesiSettori(t, k) = VIn;
        MutualInfoMetricsPaesiSettori(t, k) = MIn;
    end

    % Calculate metrics for standard modularity against partitions
    [VIn, MIn] = partition_distance(Sector, Partition{t});
    VariationInfoMetricsSettori(t, end) = VIn;
    MutualInfoMetricsSettori(t, end) = MIn;
    [VIn, MIn] = partition_distance(CountrySector, Partition{t});
    VariationInfoMetricsPaesiSettori(t, end) = VIn;
    MutualInfoMetricsPaesiSettori(t, end) = MIn;
end

%% Generate Legends for Plotting
legString = strings(length(gamma) + 1, 1);
for i = 1:length(legString)
    if i == length(legString)
        legString(i) = "Standard Modul.";
    else
        legString(i) = strcat("Modul. \gamma= ", num2str(gamma(i)));
    end
end

%% Plotting: Community Stability and Modularity Over Time

% Load aggregated community data
load CommunityAgg.mat

% Initialize arrays for storing Jaccard index and community sizes
jac = zeros(1, 15);
comSize = zeros(1, 15);
comSize2 = zeros(1, 15);

for h = 2:16
    % Process partitions to calculate community overlap and size
    pp = PartitioniOT1(:, h - 1);
    pp(isnan(pp)) = [];
    [a, ~, b] = unique(pp);
    freque = histc(b, 1:numel(a));
    [bigcomm, pos] = sort(freque);
    cc = [find(pp == a(pos(end))); find(pp == a(pos(end-1))); ...
          find(pp == a(pos(end-2))); find(pp == a(pos(end-3)))];

    pp2 = PartitioniOT1(:, h);
    pp2(isnan(pp2)) = [];
    [a2, ~, b2] = unique(pp2);
    freque = histc(b2, 1:numel(a2));
    [bigcomm2, pos2] = sort(freque);
    cc2 = [find(pp2 == a2(pos2(end))); find(pp2 == a2(pos2(end-1))); ...
           find(pp2 == a2(pos2(end-2))); find(pp2 == a2(pos2(end-3)))];
    comSize(h - 1) = length(cc2) / length(pp2);
    jac(h - 1) = length(intersect(cc, cc2)) / length(union(cc, cc2));

    pp3 = CiOT;
    [a3, ~, b3] = unique(CiOT);
    freque2 = histc(b3, 1:numel(a3));
    [bigcomm3, pos3] = sort(freque2);
    cc3 = [find(pp3 == a3(pos3(end))); find(pp3 == a3(pos3(end-1))); ...
           find(pp3 == a3(pos3(end-2))); find(pp3 == a3(pos3(end-3)))];
    comSize2(h - 1) = length(intersect(cc3, cc)) / length(union(cc3, cc));
end

% Plot community size and Jaccard index over time
figure;
yyaxis left;
plot(Tempi(2:end), comSize * 100, 'Marker', 'o', 'LineWidth', 1.5);
ylabel('Biggest 4 community size [%]');
axis tight;
grid on;

yyaxis right;
plot(Tempi(2:end), jac * 100, 'Marker', 'o', 'LineWidth', 1.5);
hold on;
plot(Tempi(2:end), comSize2 * 100, 'Marker', 'o', 'LineWidth', 1.5);
axis tight;
grid on;
ylabel('Jaccard Index [%]');
title('Community Stability');
colororder("earth");
xline(Tempi(1), '--', 'Phase 1', 'LabelVerticalAlignment', 'middle');
xline(Tempi(3), '--', 'Phase 2', 'LabelVerticalAlignment', 'middle');
xline(Tempi(8), '--', 'Phase 3', 'LabelVerticalAlignment', 'middle');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12);
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
legend('Community Size', 'Jaccard Index', 'Jaccard Index Agg.', 'Location', 'best');

%% Plot Modularity Metrics Over Time
figure;
plot(Tempi, Modul, 'Marker', 'o', 'LineWidth', 1.5);
hold on;
xline(Tempi(1), '--', 'Phase 1', 'LabelVerticalAlignment', 'middle');
xline(Tempi(4), '--', 'Phase 2', 'LabelVerticalAlignment', 'middle');
xline(Tempi(9), '--', 'Phase 3', 'LabelVerticalAlignment', 'middle');
colororder("earth");
axis tight;
grid on;
ylabel('Modularity');
title('Modularity Metric');
legend(legString, 'Location', 'best');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12);

%% Plot Mutual Information for Partitions

figure;
plot(Tempi, MutualInfoMetrics, 'Marker', 'o', 'LineWidth', 1.5);
title('Partitions');
axis tight;
grid on;
hold on;
xline(Tempi(1), '--', 'Phase 1', 'LabelVerticalAlignment', 'middle');
xline(Tempi(4), '--', 'Phase 2', 'LabelVerticalAlignment', 'middle');
xline(Tempi(9), '--', 'Phase 3', 'LabelVerticalAlignment', 'middle');
h = gca;
h.XAxis.Visible = 'off';
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12);
ylabel('Mutual Info.');
legend(legString(1:end-1), 'Location', 'best');
colororder("earth");

%% Plot Mutual Information for Country Partitions

figure;
plot(Tempi, MutualInfoMetricsPaesi, 'Marker', 'o', 'LineWidth', 1.5);
title('Countries');
axis tight;
grid on;
hold on;
h = gca;
h.XAxis.Visible = 'off';
xline(Tempi(1), '--', 'Phase 1', 'LabelVerticalAlignment', 'middle');
xline(Tempi(4), '--', 'Phase 2', 'LabelVerticalAlignment', 'middle');
xline(Tempi(9), '--', 'Phase 3', 'LabelVerticalAlignment', 'middle');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12);
colororder("earth");
ylabel('Mutual Info.');
legend(legString(1:end), 'Location', 'best');
colororder("earth");

%% Plot Mutual Information for Sector Partitions

figure;
plot(Tempi, MutualInfoMetricsSettori, 'Marker', 'o', 'LineWidth', 1.5);
title('Sectors');
axis tight;
hold on;
xline(Tempi(1), '--', 'Phase 1', 'LabelVerticalAlignment', 'middle');
xline(Tempi(4), '--', 'Phase 2', 'LabelVerticalAlignment', 'middle');
xline(Tempi(9), '--', 'Phase 3', 'LabelVerticalAlignment', 'middle');
h = gca;
h.XAxis.Visible = 'off';
grid on;
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12);
colororder("earth");
ylabel('Mutual Info.');
legend(legString(1:end), 'Location', 'best');
colororder("earth");

%% Plot Mutual Information for Country-Sector Partitions

figure;
tiledlayout(2, 1, 'TileSpacing', 'none');
nexttile;
plot(Tempi, MutualInfoMetricsPaesiSettori, 'Marker', 'o', 'LineWidth', 1.5);
title('Country-Sectors');
axis tight;
hold on;
xline(Tempi(1), '--', 'Phase 1', 'LabelVerticalAlignment', 'middle');
xline(Tempi(4), '--', 'Phase 2', 'LabelVerticalAlignment', 'middle');
xline(Tempi(9), '--', 'Phase 3', 'LabelVerticalAlignment', 'middle');
h = gca;
h.XAxis.Visible = 'off';
grid on;
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12);
colororder("earth");
ylabel('Mutual Info.');
nexttile;
plot(Tempi, VariationInfoMetricsPaesiSettori, 'Marker', 'o', 'LineWidth', 1.5);
axis tight;
grid on;
hold on;
xline(Tempi(1), '--', 'Phase 1', 'LabelVerticalAlignment', 'middle');
xline(Tempi(4), '--', 'Phase 2', 'LabelVerticalAlignment', 'middle');
xline(Tempi(9), '--', 'Phase 3', 'LabelVerticalAlignment', 'middle');
ylabel('Variation Info.');
legend(legString(1:end), 'Location', 'best');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12);
colororder("earth");

%% Display Partition Information

% Calculate the number of communities in each partition
numcomm = cellfun(@max, Partition);
numcommOT = cellfun(@max, PartitionOT);
Colord = orderedcolors("earth");

figure
tiledlayout(4, 1, 'TileSpacing', 'none');
nexttile;
violin(Partition', 'facealpha', 1, 'facecolor', Colord(4,:), 'edgecolor', Colord(4,:));
axis tight;
grid on;
xticks(1:length(Tempi));
xticklabels(strings(size(Tempi)));
ylabel('Num. C.');
TextLocation('Standard Modularity', 'Location', 'best');

nexttile;
violin(PartitionOT(:, 1)', 'facealpha', 1, 'facecolor', Colord(1,:), 'edgecolor', Colord(1,:), 'bw', .5);
axis tight;
grid on;
ylabel('Num. C.');
TextLocation('Modularity \gamma=0.0001', 'Location', 'best');
xticks(1:length(Tempi));
xticklabels(strings(size(Tempi)));

nexttile;
violin(PartitionOT(:, 2)', 'facealpha', 1, 'facecolor', Colord(2,:), 'edgecolor', Colord(2,:));
axis tight;
grid on;
ylabel('Num. C.');
TextLocation('Modularity \gamma=0.001', 'Location', 'best');
xticks(1:length(Tempi));
xticklabels(strings(size(Tempi)));

nexttile;
violin(PartitionOT(:, 3)', 'facealpha', .4, 'facecolor', Colord(3,:), 'edgecolor', Colord(3,:));
axis tight;
grid on;
ylabel('Num. C.');
TextLocation('Modularity \gamma=0.01', 'Location', 'best');
xticks(1:length(Tempi));
xticklabels(Tempi);

sgtitle('Community Distributions');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18);
