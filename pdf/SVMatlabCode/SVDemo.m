
%% Demo: Construct the spot volatility index (SV) 
% Note: make sure the sample data DemoData.mat and function GetV.m are in
% the current working path

clc;clear;
% load the sample data
load('DemoData.mat');

% call GetV to calculate the option-implied spot variance 
% based on the near and next term options in Equation (2)
Vnear = GetV(options1,T1,spot);
Vnext = GetV(options2,T2,spot);

% average and rescale to get the SV Index in Equation (1)
SV = 100*sqrt((Vnear + Vnext)/2);