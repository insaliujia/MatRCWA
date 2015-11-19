% Rigorous Coupled-Wave Analysis
% The program is using object oriented program at MATLAB
% Jia LIU Ph.D student, INSA de Lyon
% Supervisor Regis Orobtchouk, INSA de Lyon
% This RCWA method is modified to modle silicon photovoltaics devices 
%% INITIALIZE MATLAB 
close all; 
clear all; 
clear classes;
addpath(genpath('C:\Users\jliu\RCWA'));
clc;
%%
%% Build new simulation object, source and device
ShowProcess = 1;    % if you cannot delete it try this: set(0,'ShowHiddenHandles','on'); delete(get(0,'Children'))
% Build new simulation input: ref index [er,ur], trn index [er,ur],waitbar
Simul = RCWA([1,1],[1,1],ShowProcess);
% Build source input: [wavelength],[theta(angle with axis XZ),phi(angle in
% the XY plane)], polarization
S = Source([300:1100],[0,0],[1/sqrt(2),1i/sqrt(2)]);
Air = Material('Air',[1,1]);
Si = Material('silicon');
Si.Selectn(S);
% Build device input: [length,width],[number in x,number in y](should be in proportion),[spatial harmonics]
D = Device([0.3,0.3],[1500,1],[65,1]);
% AddLayer(D,Si,0.15,1);
AddLayer(D,Si,0.1,1);
AddLayer(D,Si,0.9,1);
AddPattern(D,'Rectangle',[0.3/2,0.3/2],[0.3*0.5,0.3],1,Air);
% AddPattern(D,'Rectangle',[0.15,0.2],[0.1,0.4],1,Air);
% AddPattern(D,'Rectangle',[0.35,0.2],[0.1,0.4],1,Air);
% Simulation
RCWARun(Simul,S,D);
PlotRT(Simul);