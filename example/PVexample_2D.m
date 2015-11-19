% Rigorous Coupled-Wave Analysis
% The program is using object oriented program at MATLAB
% Jia LIU Ph.D student, INSA de Lyon
% Supervisor Regis Orobtchouk, INSA de Lyon
% This RCWA method is modified to modle silicon photovoltaics devices 
%% INITIALIZE MATLAB 
close all; 
clear all; 
clear classes;
addpath(genpath('C:\Users\jliu\RCWA'))
clc;
%%
%% Build new simulation object, source and device
ShowProcess = 1;    % if you cannot delete it try this: set(0,'ShowHiddenHandles','on'); delete(get(0,'Children'))
% Build new simulation input: ref index [er,ur], trn index [er,ur],waitbar
Simul = RCWA([1,1],[1,1],ShowProcess);
% UseShrinkEffect(Simul,4);
% Record diffrection order
% RecordDiffOrder(Simul);
% Build source input: [wavelength],[theta(angle with axis XZ),phi(angle in the XY plane)], polarization
S = Source([400:410],[0,0],[1/sqrt(2),1/sqrt(2)]);
% Ag = Material('Ag');
Si = Material('silicon');
Air = Material('Air',1);
% Al = Material('Al');
Si.Selectn(S);
% Build device input: [length,width],[number in x,number in y](should be in proportion),[spatial harmonics]
D = Device([0.06,0.06],[1024,1024],[25,25]);
% AddLayer(D,Air,0.1,1);
AddLayer(D,Si,0.06,1);
% AddPattern(D,'Cylinder',[0.3,0.3],0.2,1,Al);
% AddPattern(D,'Cylinder',[0.3,0.3],0.15,1,Air);
% Dispersion(Simul,0);
% D.BuildLayer(1);
% D.BuildPattern(1);
% AddPattern(D,'Cylinder',[0.5,0.5],0.15,[1:2],Al);
% AddPattern(D,'Rectangle',[0.5,0.5],[0.1,0.3],[1:4],Si)
% AddPattern(D,'Rectangle',[0.3,0.6],[0.1,0.2],[1:4],Si);
% AddPattern(D,'Triangle',[0.3,0.6],0.15,[1:5],ITO);
% AddPattern(D,'Cylinder',[0.45,0.15],0.09,1,Si);
AddPattern(D,'Cylinder',[0.03,0.03],0.02,1,Air);
% AddPattern(D,'Cylinder',[0.15,0.15],0.05,1,Si);
% AddPattern(D,'Cylinder',[0.45,0.45],0.05,1,Si);
% AddPattern(D,'Rectangle',[0.3,0.3],[0.1,0.2],4,Air);
% AddPattern(D,'Pyramid',[0.3,0.3],[0.3,0.3],[1:5],Air);
RCWARun(Simul,S,D)
PlotRT(Simul)


