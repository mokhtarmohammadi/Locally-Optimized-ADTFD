%%  Author:
%     Mokhtar Mohammadi
%In this code we assume that the user add in this path the TFSAP toolboox.
% Please cite "Mokhtar Mohammadi, Ali Akbar Pouyan, Nabeel Ali Khan, Vahid Abolghasemi:
%              Locally Optimized Adaptive Directional Time-Frequency Distributions. 
%               CSSP 37(8): 3154-3174 (2018)" 

%% Locating and Adding Functions Directory
clear all;
close all;
currentDirectory = pwd;
[upperPath, ~, ~] = fileparts(currentDirectory);
addpath([upperPath '\Core'])
addpath([upperPath '\TFSA7'])
addpath([upperPath '\Core\AOK'])
addpath([upperPath '\TFTB2'])

%% Signal Parameters and Loading
n = 0:255; fs = 32; tr = 3;
load sigEEGseizure;

%% TFD Generation
TFD_WVD = quadtfd( sig, length(sig)-1, 1, 'wvd');
TFD_CKD  = cmpt(sig, 'ckd', 1, 0.075, 0.075);
TFD_EMBD = quadtfd(sig, length(sig)-1, 1, 'emb', 0.1, 0.2);
% TFD_SM   = specSM(sig, 4, 85, 'hann', 84, 512);
TFD_SM=Adaptive_S_method(sig,21,'hann');
 TFD_RGK = rgk(sig,2);
[TFD_FADTFD orient]=FADTFD_EEGR(sig);
TFD_SPEC = quadtfd( sig, length(sig)-1, 1, 'specx', 21, 'hann',256);
I_max_new=adaptive_optimal_tfd(sig);
TFD_AOK = real(I_max_new);%% Plotting
figure; SetFigDef(16,9,'Times',20); tfsapl(sig,TFD_WVD(:,1:tr:end),'SampleFreq',fs,'Title','(a)','TimePlot', 'on', 'FreqPlot', 'on', 'TFfontSize' , 20,'plotfn','pcolor'); 
figure; SetFigDef(16,9,'Times',20);tfsapl(sig,TFD_SPEC(:,1:tr:end),'SampleFreq',fs,'Title','(b)','TimePlot', 'on', 'FreqPlot', 'on', 'TFfontSize' , 20,'plotfn','pcolor');
figure; SetFigDef(16,9,'Times',20);tfsapl(sig,TFD_EMBD(:,1:tr:end),'SampleFreq',fs,'Title','(c)','TimePlot', 'on', 'FreqPlot', 'on', 'TFfontSize' , 20,'plotfn','pcolor');
figure; SetFigDef(16,9,'Times',20);tfsapl(sig,TFD_SM(:,1:tr:end),'SampleFreq',fs,'Title','(d)','TimePlot', 'on', 'FreqPlot', 'on', 'TFfontSize' , 20,'plotfn','pcolor');
figure; SetFigDef(16,9,'Times',20);tfsapl(sig,TFD_CKD(:,1:tr:end),'SampleFreq',fs,'Title','(e)','TimePlot', 'on', 'FreqPlot', 'on', 'TFfontSize' , 20,'plotfn','pcolor');
figure; SetFigDef(16,9,'Times',20);tfsapl(sig,TFD_AOK(:,1:tr:end),'SampleFreq',fs,'Title','(f)','TimePlot', 'on', 'FreqPlot', 'on', 'TFfontSize' , 20,'plotfn','pcolor');
figure; SetFigDef(16,9,'Times',20);tfsapl(sig,TFD_RGK(:,1:tr:end),'SampleFreq',fs,'Title','(g)','TimePlot', 'on', 'FreqPlot', 'on', 'TFfontSize' , 20,'plotfn','pcolor');
figure; SetFigDef(16,9,'Times',20);tfsapl(sig,TFD_FADTFD(:,1:tr:end),'SampleFreq',fs,'Title','(h)','TimePlot', 'on', 'FreqPlot', 'on', 'TFfontSize' , 20,'plotfn','pcolor');
% t=0:1/32:8-1/32;
% figure;plot(t,sig);