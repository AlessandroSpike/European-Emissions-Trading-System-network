%% DistanceAnalysisReal
% This script analyzes trade/interaction patterns based on geographical distance
% and country borders using real data. It processes trade amounts between
% different locations and analyzes them across distance quartiles and borders.
%
% Input required:
% - OutputGravity.mat: Contains trade data (iso_d, iso_o, amount, year)
% - AdjMatrixTime.mat: Contains time periods, location codes, and distance matrix
% - Border_Distance_Info: Contains border information and location names
%
% Output:
% - ConfiniReal.mat: Contains analyzed trade patterns by distance quartiles
%   and border regions

clc;clearvars;close all

%% Add data path
addpath("Initialdata")

%% Load and process input data
% Load trade/interaction data
load('OutputGravity.mat','iso_d','iso_o',...
    'amount','year')
iso_d=string(iso_d.iso_d(2:end));      % Destination locations
iso_o=string(iso_o.iso_o(2:end));      % Origin locations
amount=amount.amount(2:end);            % Trade/interaction amounts
year=year.year(2:end);                 % Time periods

% Load geographical and border information
load('AdjMatrixTime.mat','Tempi','C','Distance')
load('Border_Distance_Info','Nomi','DummyBorder');
NomiBorder=Nomi;
Call=C;

%% Create time-series adjacency matrices
AdjT=zeros(length(Call),length(Call),length(Tempi));
for t=1:length(Tempi)
    % Extract data for current time period
    tengo=year==(Tempi(t));
    
    % Process amounts (replace -1 with small value)
    amount_t=abs(amount(tengo));
    amount_t(amount_t==-1)=0.001;
    iso_ot=iso_o(tengo);
    isodt=iso_d(tengo);
    
    % Create adjacency matrix for current period
    Nomi=[iso_ot,isodt];
    [C,~,ib] = unique(Nomi, 'stable');
    ib = reshape(ib, size(Nomi));
    Adj=full(sparse(ib(:,1),ib(:,2),amount_t,length(C),length(C)));
    
    % Expand matrix to full size
    n_large = length(Call);
    large_matrix = zeros(n_large);
    [~, small_indices] = ismember(C, Call);
    [I, J] = meshgrid(small_indices, small_indices);
    large_matrix(sub2ind(size(large_matrix), I(:), J(:))) = Adj(:);
    
    AdjT(:,:,t)=large_matrix;
end

%% Analyze by distance quartiles
Amount=sum(AdjT,3);
Amount=Amount(:);
Distanza=Distance(:);

% Create pairs of locations
n = size(large_matrix, 1);   
[row_indices, col_indices] = meshgrid(1:n, 1:n);
row_indices = row_indices(:);
col_indices = col_indices(:);  
Nomi = [Call(row_indices), Call(col_indices)];

% Split into distance quartiles
dist1=find(Distanza<prctile(Distanza,25));
dist2=find(Distanza>prctile(Distanza,25) & Distanza<prctile(Distanza,50));
dist3=find(Distanza>prctile(Distanza,50) & Distanza<prctile(Distanza,75));
dist4=find(Distanza>prctile(Distanza,75));

% Extract amounts for each distance quartile
Amount1=Amount(dist1);
Amount2=Amount(dist2);
Amount3=Amount(dist3);
Amount4=Amount(dist4);

% Analyze domestic vs international patterns
Nomi2=extractBefore(Nomi,3);
StessoPaese=Nomi2(:,1)==Nomi2(:,2);
StessoPaese1=StessoPaese(dist1);
StessoPaese2=StessoPaese(dist2);
StessoPaese3=StessoPaese(dist3);
StessoPaese4=StessoPaese(dist4);

% Split amounts by domestic/international for each quartile
Amount1a=Amount1(StessoPaese1==0);  % International
Amount1b=Amount1(StessoPaese1==1);  % Domestic
Amount2a=Amount2(StessoPaese2==0);
Amount2b=Amount2(StessoPaese2==1);
Amount3a=Amount3(StessoPaese3==0);
Amount3b=Amount3(StessoPaese3==1);
Amount4a=Amount4(StessoPaese4==0);
Amount4b=Amount4(StessoPaese4==1);

%% Analyze border regions
% Expand border dummy variable to full size
DummyBorder2 = zeros(size(Nomi,1),1);    
[~, indices] = ismember(join(NomiBorder), join(Nomi));
DummyBorder2(indices) = DummyBorder;

% Analyze trade patterns in border regions
TileConfine=Nomi(DummyBorder2==1,:);
TileConfine=unique(TileConfine(:));

% Calculate domestic and international fractions for border regions
frazione=zeros(size(TileConfine));    % For domestic trade
frazione2=zeros(size(TileConfine));   % For international trade
for t=1:length(TileConfine)
    t
    chi=TileConfine(t);
    chiS=extractBefore(chi,3);
    chi1=find(Nomi(:,1)==chi);
    chi1b=Nomi2(chi1,2);
    chi2=find(Nomi(:,2)==chi);
    chi2b=Nomi2(chi2,1);

    % Calculate domestic and international trade volumes
    uno=Amount(chi1);
    due=Amount(chi2);
    
    uno_s=sum(uno(chi1b==chiS));
    uno_d=sum(uno(chi1b~=chiS));
    due_s=sum(due(chi2b==chiS));
    due_d=sum(due(chi2b~=chiS));

    frazione(t)=(uno_s+due_s);      % Domestic trade
    frazione2(t)=(uno_d+due_d);     % International trade
end

% Clean and process border region data
frazione(isinf(frazione))=nan;
frazione2(isinf(frazione2))=nan;
TileConfine2=extractBefore(TileConfine,3);
[G,paesi] = findgroups(TileConfine2);
varBP = splitapply(@nanmean,frazione,G);
varBP2 = splitapply(@nanmean,frazione2,G);
paesi(varBP==0)=[];
varBP(varBP==0)=[];
varBP2(varBP2==0)=[];

%% Save results
save('ConfiniReal.mat','Amount1a','Amount2a','Amount3a','Amount4a',...
    'Amount1b','Amount2b','Amount3b','Amount4b','varBP','varBP2','paesi')