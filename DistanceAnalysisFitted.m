%% DistanceAnalysisFitted
% This script analyzes trade/interaction patterns based on geographical distance
% and country borders using fitted values from a gravity model. It processes 
% fitted interaction values between different locations and analyzes them 
% across distance quartiles and borders.
%% Initialize Environment
% Clear workspace and figures for clean execution
clc;clearvars;close all

%% Set Data Path
% Add path containing input data files
addpath("Initialdata")
%% Load Input Data
% Load gravity model outputs including ISO codes and fitted values
% OutputGravity.mat contains:
% - iso_d: Destination ISO codes
% - iso_o: Origin ISO codes
% - fitted: Fitted values from gravity model
% - year: Time periods
load('OutputGravity.mat','iso_d','iso_o','fitted','year')
iso_d=string(iso_d.iso_d(2:end));  % Convert destination ISO codes to string
iso_o=string(iso_o.iso_o(2:end));  % Convert origin ISO codes to string
fitted=fitted.xfittedvalues(2:end); % Get fitted values
year=year.year(2:end);             % Get years

% Load temporal adjacency matrix and geographical information
% AdjMatrixTime.mat contains:
% - Tempi: Time periods
% - C: Location codes
% - Distance: Distance matrix between locations
load('AdjMatrixTime.mat','Tempi','C','Distance')

% Load border information
% Border_Distance_Info contains:
% - Nomi: Location names
% - DummyBorder: Binary indicator for border regions
load('Border_Distance_Info','Nomi','DummyBorder');
NomiBorder=Nomi;
Call=C;

%% Construct Temporal Adjacency Matrices
% Create 3D matrix of interactions over time (nodes x nodes x time)
AdjT=zeros(length(Call),length(Call),length(Tempi));

% Process each time period
for t=1:length(Tempi)
    % Filter data for current time point
    tengo=year==(Tempi(t));
    
    % Extract relevant data for current time
    fitted_t=fitted(tengo);      % Fitted values
    iso_ot=iso_o(tengo);        % Origin ISO codes
    isodt=iso_d(tengo);         % Destination ISO codes
    Nomi=[iso_ot,isodt];        % Combined location pairs
    
    % Create unique node list and get indices
    [C,~,ib] = unique(Nomi, 'stable');
    ib = reshape(ib, size(Nomi));
    
    % Create adjacency matrix for current time point
    Adj=full(sparse(ib(:,1),ib(:,2),fitted_t,length(C),length(C)));
    
    % Map to full network size
    n_large = length(Call);
    large_matrix = zeros(n_large);
    [~, small_indices] = ismember(C, Call);
    [I, J] = meshgrid(small_indices, small_indices);
    large_matrix(sub2ind(size(large_matrix), I(:), J(:))) = Adj(:);
    
    % Store in temporal adjacency tensor
    AdjT(:,:,t)=large_matrix;
end

%% Distance-based Analysis
% Calculate total interaction amounts and prepare distance data
Amount=sum(AdjT,3);
Amount=Amount(:);
Distanza=Distance(:);

% Create node pairs for analysis
n = size(large_matrix, 1);   
[row_indices, col_indices] = meshgrid(1:n, 1:n);
row_indices = row_indices(:);
col_indices = col_indices(:);  
Nomi = [Call(row_indices), Call(col_indices)];

% Divide distances into quartiles for analysis
dist1=find(Distanza<prctile(Distanza,25));                  % 0-25th percentile
dist2=find(Distanza>prctile(Distanza,25) & Distanza<prctile(Distanza,50));  % 25-50th percentile
dist3=find(Distanza>prctile(Distanza,50) & Distanza<prctile(Distanza,75));  % 50-75th percentile
dist4=find(Distanza>prctile(Distanza,75));                  % 75-100th percentile

% Group interaction amounts by distance quartiles
Amount1=Amount(dist1);  % Shortest distance
Amount2=Amount(dist2);  % Short-medium distance
Amount3=Amount(dist3);  % Medium-long distance
Amount4=Amount(dist4);  % Longest distance

%% Same-Country Analysis
% Identify pairs within the same country using ISO code prefixes
Nomi2=extractBefore(Nomi,3);
StessoPaese=Nomi2(:,1)==Nomi2(:,2);  % Binary indicator for same country

% Split same-country indicator by distance quartiles
StessoPaese1=StessoPaese(dist1);
StessoPaese2=StessoPaese(dist2);
StessoPaese3=StessoPaese(dist3);
StessoPaese4=StessoPaese(dist4);

% Separate amounts by distance quartile and country relationship
% Different countries (suffix 'a')
Amount1a=Amount1(StessoPaese1==0);
Amount2a=Amount2(StessoPaese2==0);
Amount3a=Amount3(StessoPaese3==0);
Amount4a=Amount4(StessoPaese4==0);

% Same country (suffix 'b')
Amount1b=Amount1(StessoPaese1==1);
Amount2b=Amount2(StessoPaese2==1);
Amount3b=Amount3(StessoPaese3==1);
Amount4b=Amount4(StessoPaese4==1);

%% Border Region Analysis
% Expand border dummy variable to full pair matrix
DummyBorder2 = zeros(size(Nomi,1),1);    
[~, indices] = ismember(join(NomiBorder), join(Nomi));
DummyBorder2(indices) = DummyBorder;

% Identify border regions
TileConfine=Nomi(DummyBorder2==1,:);
TileConfine=unique(TileConfine(:));

% Initialize arrays for within-country and cross-border interactions
frazione=zeros(size(TileConfine));    % Within-country interactions
frazione2=zeros(size(TileConfine));   % Cross-border interactions

% Calculate interactions for each border region
for t=1:length(TileConfine)
    t  % Display progress
    chi=TileConfine(t);               % Current region
    chiS=extractBefore(chi,3);        % Country code
    
    % Find connections where current region is origin/destination
    chi1=find(Nomi(:,1)==chi);        % Region as origin
    chi1b=Nomi2(chi1,2);             % Destination countries
    chi2=find(Nomi(:,2)==chi);        % Region as destination
    chi2b=Nomi2(chi2,1);             % Origin countries

    % Get interaction amounts
    uno=Amount(chi1);                 % Outgoing interactions
    due=Amount(chi2);                 % Incoming interactions
    
    % Sum domestic and international interactions
    uno_s=sum(uno(chi1b==chiS));      % Domestic outgoing
    uno_d=sum(uno(chi1b~=chiS));      % International outgoing
    due_s=sum(due(chi2b==chiS));      % Domestic incoming
    due_d=sum(due(chi2b~=chiS));      % International incoming

    % Store total domestic and international interactions
    frazione(t)=(uno_s+due_s);        % Total domestic
    frazione2(t)=(uno_d+due_d);       % Total international
end

% Clean infinite values
frazione(isinf(frazione))=nan;
frazione2(isinf(frazione2))=nan;

% Calculate country-level averages
TileConfine2=extractBefore(TileConfine,3);
[G,paesi] = findgroups(TileConfine2);
varBP = splitapply(@nanmean,frazione,G);     % Average domestic interactions
varBP2 = splitapply(@nanmean,frazione2,G);   % Average international interactions

% Remove countries with no interactions
tolgo=sum([varBP, varBP2],2)==0;
paesi(tolgo)=[];
varBP(tolgo)=[];
varBP2(tolgo)=[];
paesiFitted=paesi;

%% Save Results
% Prepare fitted values for export
Amount1afitted=Amount1a;  % Cross-country, shortest distance
Amount2afitted=Amount2a;  % Cross-country, short-medium distance
Amount3afitted=Amount3a;  % Cross-country, medium-long distance
Amount4afitted=Amount4a;  % Cross-country, longest distance

Amount1bfitted=Amount1b;  % Same-country, shortest distance
Amount2bfitted=Amount2b;  % Same-country, short-medium distance
Amount3bfitted=Amount3b;  % Same-country, medium-long distance
Amount4bfitted=Amount4b;  % Same-country, longest distance

varBPfitted=varBP;       % Average domestic border region interactions
varBP2fitted=varBP2;     % Average international border region interactions

% Save analyzed data to mat file
save('ConfiniFitted.mat','Amount1afitted','Amount2afitted','Amount3afitted',...
    'Amount4afitted','Amount1bfitted','Amount2bfitted','Amount3bfitted',...
    'Amount4bfitted','varBPfitted','varBP2fitted','paesiFitted')