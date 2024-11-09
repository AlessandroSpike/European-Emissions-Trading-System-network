% Description:
% This script visualizes the European Union Emissions Trading System (EU ETS) network
% using geographical data and adjacency matrices.
% 
% Input:
%   - NUTS_RG_60M_2021_3035.shp: Shape file containing EU regional boundaries
%   - AdjMatrixAggregate.mat: Adjacency matrix representing connections between regions
%
% Output:
%   - Geographical visualization of the EU ETS network showing:
%     * Nodes: Regions sized by their connection strength
%     * Edges: Connections between regions
%     * Colors: Based on outgoing connection strength (logarithmic scale)

% Clear workspace and figures
clc;clearvars;close all

%% Add required paths
addpath('NUTS_RG_60M_2021_3035A.shp')
addpath('Initialdata')

%% Load adjacency matrix data
load AdjMatrixAggregate.mat

%% Load and process shape file data
S = shaperead('NUTS_RG_60M_2021_3035.shp');
Slabel=string({S.NUTS_ID});

% Update region labels to match current NUTS classification
% These replacements handle changes in NUTS region codes
Slabel=replace(Slabel,'BE224','BE221');
% ... [remaining replacements] ...

% Find matching regions between shape file and adjacency matrix
[chi,tengo,ia]=intersect(Slabel,C);
S=S(tengo);
AdjN=AdjN(ia,:);
AdjN=AdjN(:,ia);
C=C(ia);

% Extract and calculate centroid coordinates
X={S.X};
Y={S.Y};
x=cellfun(@nanmean,X);
y=cellfun(@nanmean,Y);

% Convert projected coordinates to latitude/longitude
info = shapeinfo('NUTS_RG_60M_2021_3035.shp');
proj = info.CoordinateReferenceSystem;
[Latit,Longit] = projinv(proj,x,y);

% Remove regions outside the area of interest
tolgo1=Latit<35 | Latit>70;  % Latitude bounds
tolgo2=Longit<-30 | Longit>30;  % Longitude bounds
tolgo3=tolgo1+tolgo2;
tolgo=find(tolgo3>0);
Latit(tolgo)=[];
Longit(tolgo)=[];
AdjN(tolgo,:)=[];
AdjN(:,tolgo)=[];
C(tolgo)=[];

% Calculate in-degree and out-degree strengths
inS=sum(AdjN);
inS(inS==0)=.1;  % Avoid zero values
outS=sum(AdjN,2);
outS(outS==0)=.1;  % Avoid zero values

% Sort regions by out-degree strength
[outS,pos]=sort(outS);
AdjN=AdjN(pos,:);
AdjN=AdjN(:,pos);
C=C(pos);
Latit=Latit(pos);
Longit=Longit(pos);

% Generate color map based on logarithmic out-degree strength
rgb = vals2colormap(log(outS), 'abyss');

% Prepare edge weights for plotting
pesi = reshape(AdjN,size(AdjN,1)^2,1);
pesi(pesi==0)=[];
soglia = 0;  % Threshold for displaying connections

%% Create visualization
figure
% Plot edges (connections between regions)
for x = 1:length(AdjN)
    for y = 1:length(AdjN)
        l = AdjN(x,y);
        if l > soglia
            p1=extractBefore(C(x),3);  % Get country code
            p2=extractBefore(C(y),3);
            if sum(isnan([Longit(y) Latit(y) Longit(x) Latit(x)]))==0
                if p1~=p2
                    % Inter-country connections in black
                    geoplot([Latit(y) Latit(x)],[Longit(y) Longit(x)],'k:','LineWidth',0.05)
                else
                    % Intra-country connections in red
                    geoplot([Latit(y) Latit(x)],[Longit(y) Longit(x)],'r:','LineWidth',0.05)
                end
                hold on
            end
        end
    end
end

% Plot nodes (regions)
for x = 1:length(AdjN)
    for y = 1:length(AdjN)
        l = AdjN(x,y);
        if l > soglia
            if sum(isnan([Longit(y) Latit(y) Longit(x) Latit(x)]))==0
                % Node size based on in-degree strength
                geoscatter(Latit(y),Longit(y),20*(inS(y)/max(inS)+.1),rgb(y,:),'filled')
                geoscatter(Latit(x),Longit(x),20*(inS(x)/max(inS))+.1,rgb(x,:),'filled')
                hold on
            end
        end
    end
end

% Set map properties
geolimits([35 70],[-15 30])  % Set geographical bounds
title('Aggregate EU ETS network')
geobasemap darkwater  % Set base map style