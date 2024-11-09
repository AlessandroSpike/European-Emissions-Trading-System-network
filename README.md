# European-Emissions-Trading-System-network
This folders contains files to reproduce results of the paper "Carbon trade biases and the emerging mesoscale structure of the European Emissions Trading System network"
by Andrea Flori and Alessandro Spelta.

%% GRAVITY ESTIMATION
File name: gravity_model.R
This file reproduces the results reported in Table 1 and in Supplementary Information S1. The code is based on R.
The analysis requires a set of R packages to be installed. They are reported in the first 38 lines.

Line 46: we import the main .csv file, where relevant information is located.
Lines 49-51: we create the Same Country (C) variable
Lines 54-69: we create the Same Sector (S) variable for each year (as a function of aggregate verified emissions of the installations in the sample in each NUTS-3 province)
Lines 73-215: for each year separately, the code computes the gravity model estimation using the package ppml. 
Specifically, the y variable is the amount of traded allowances in each year. The geodesic distance is labelled "distance". The pedix "acq" refers to the destination node, while the pedix "transf" to the origin node. We fit Eq. 1 of the paper.

Lines 219-509: we construct summary statistics for the estimates of the gravity models. These estimates are conveniently saved in the object "stats_db".
Lines 513-535: we create the database for the plots.
Lines 538-735: these lines build and export all the plots reported in the paper for the main analysis.

Lines 740-790: we calculate the deviance and the fitted values for each model. 

%% COMMUNITY DETECTION ANALYSIS
%BiasCountrySectorReal
This script calculates three types of bias (country, sector, country-sector)
in both weighted and binary network connections over time. Bias is measured
as percentage of within-group connections relative to total connections.
Input required: AdjMatrixTime.mat containing adjacency matrices (Adj) and 
timestamps (Tempi), plus partition vectors for countries and sectors.
Output: BiasReal.mat containing bias calculations over time.
%BiasCountrySectorFitted
This script analyzes network bias in gravity model outputs across three dimensions:
country, sector, and country-sector combinations. For each dimension, both weighted
and binary (unweighted) bias metrics are calculated.
Input required: 
- OutputGravity.mat: Contains gravity model results (iso_d, iso_o, fitted, year)
- AdjMatrixTime.mat: Contains partition vectors for countries and sectors
Output: 
- BiasFitted.mat: Contains bias calculations (biasFitted) for each time period (Tempi)
The bias is calculated as percentage of within-group connections relative to
total connections, for both weighted and unweighted network representations.
%DistanceAnalysisReal
This script analyzes trade/interaction patterns based on geographical distance
and country borders using real data. It processes trade amounts between
different locations and analyzes them across distance quartiles and borders.
Input required:
- OutputGravity.mat: Contains trade data (iso_d, iso_o, amount, year)
- AdjMatrixTime.mat: Contains time periods, location codes, and distance matrix
- Border_Distance_Info: Contains border information and location names
Output:
- ConfiniReal.mat: Contains analyzed trade patterns by distance quartiles
  and border regions
%DistanceAnalysisFitted
This script analyzes trade/interaction patterns based on geographical distance
and country borders using fitted values from a gravity model. It processes 
fitted interaction values between different locations and analyzes them 
across distance quartiles and borders.
Input required:
- OutputGravity.mat: Contains fitted gravity model values and ISO codes
- AdjMatrixTime.mat: Contains time periods, location codes, and distance matrix
- Border_Distance_Info: Contains border information and location names
Output:
- ConfiniFitted.mat: Contains analyzed interaction patterns by distance
  quartiles and border regions
%PlotBias
This script creates visualizations comparing real and fitted network characteristics,
including bias analysis, border effects, and distance-based weight distributions.
Input required:
- BiasReal.mat: Contains bias measures from real network
- BiasFitted.mat: Contains bias measures from fitted network
- ConfiniFitted.mat: Contains analyzed fitted patterns by distance and borders
- ConfiniReal.mat: Contains analyzed real patterns by distance and borders
Outputs:
- Figure 1: Bias comparison and border region analysis
- Figure 2: Distance-based weight distribution analysis
%GeoNetPlot
This script visualizes the European Union Emissions Trading System (EU ETS) network
using geographical data and adjacency matrices.
Input:
  - NUTS_RG_60M_2021_3035.shp: Shape file containing EU regional boundaries
  - AdjMatrixAggregate.mat: Adjacency matrix representing connections between regions
Output:
  - Geographical visualization of the EU ETS network showing:
    * Nodes: Regions sized by their connection strength
    * Edges: Connections between regions
    * Colors: Based on outgoing connection strength (logarithmic scale)
%PlotDeviance
This MATLAB script analyzes and visualizes deviance data for different gravity models
of carbon emission trading. It compares various model configurations and
calculates confidence intervals for the difference in deviances.
Input:
- CSV files containing deviance data for different model configurations
- Years range from 2005 to 2020
Output:
- Two figures: 
  1. Deviance plots for different model configurations
  2. Log of deviance difference with confidence intervals
- Console output of difference in deviances and confidence intervals
%CommunityOT
This script performs a network analysis over time using Optimal Transport
(OT) methods and classical modularity approaches. It computes modularity
and partition metrics for a set of temporal networks and compares them
with community structures based on countries and sectors.
Inputs:
- OutputGravity.mat: Contains data on countries, fitted values, and amounts.
- AdjMatrixTime.mat: Contains adjacency matrices and partition information.
- CommunityAgg.mat: Provides aggregated community information.
Outputs:
- Plots and figures representing modularity and partition metrics over time.
- Saved figures in 'PlotFig/' directory.
%CommunityOT_AggregateNet
This script analyzes transportation and cultural distance data across
different regions over time. It uses optimal transport methods to assess
community interactions based on given data. The script outputs plots
showing average cultural distances within and between communities.
Inputs:
- 'OutputGravity.mat': Contains:
    - iso_d: Destination identifiers.
    - iso_o: Origin identifiers.
    - amount: Amount of transport or interaction.
    - fitted: Fitted values from a gravity model.
    - year: Corresponding years for the data.
- 'AdjMatrixTime.mat': Contains:
    - Tempi: Time periods for the analysis.
    - C: List of community identifiers.
    - Distance: Matrix of distances (not used in this script).
- 'Border_Distance_Info.mat': Contains:
    - Nomi: Border information (not used in this script).
    - DummyBorder: Dummy variables for border information (not used).
- 'culturaldistance.mat': Contains:
    - ctr: Country codes.
    - Various columns of cultural distance data.
Outputs:
- Figures showing the average cultural distances within and between
  communities based on different factors.
- Optional: Saves community aggregation data to 'CommunityAgg.mat'.

