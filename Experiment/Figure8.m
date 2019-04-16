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
fs = 16;tr=1;
n = 0:1/fs:8-1/fs; N = length(n);
sig1 = 2*cos(2*pi*n*0.05*fs);
sig2 = 2*cos(2*pi*n*0.1*fs);
sigm = 0.0020;
B = fir1(16*2,0.08,'high');
s1 = 10*exp(-(fs^2*(n-0.9375).^2)/sigm);
s1 = filter(B,1,s1);
s2 = 10*exp(-(fs^2*(n-2.1875).^2)/sigm);%.*cos(0.1*pi*n);
s2 = filter(B,1,s2);
s3 = 10*exp(-(fs^2*(n- 3.4375).^2)/sigm);%.*cos(0.1*pi*n);
s3 = filter(B,1,s3);
s4 = 10*exp(-(fs^2*(n-4.6875).^2)/sigm);%.*cos(0.1*pi*n);
s4 = filter(B,1,s4);
s5 = 10*exp(-(fs^2*(n- 5.9375).^2)/sigm);%.*cos(0.1*pi*n);
s5 = filter(B,1,s5);
s = s1 + s2 + s3 + s4 + s5;
Signal = sig1 + sig2 + s;
sig=Signal;
Signal = awgn(Signal,10,'measured'); 
load simsig;
%% True TFD
I1 = wvd(s1, N-1) + wvd(s2, N-1) + wvd(s3, N-1) + wvd(s4, N-1) + wvd(s5, N-1);
I2 = wvd(sig1, N-1);
I3 = wvd(sig2, N-1);
TFD_True = I1 + I2 + I3;

%% TFD Generation
TFD_SPEC = quadtfd( Signal, length(Signal)-1, 1, 'specx', 23, 'hann',256);
TFD_EMBD = quadtfd(Signal,length(Signal)-1, 1, 'emb', 0.1, 0.1);
TFD_CKD = cmpt( Signal, 'ckd', 1, 0.1, 0.1);
TFD_SM   = Adaptive_S_method(Signal,23,'hann');
I_max_new=adaptive_optimal_tfd(Signal);
TFD_AOK = real(I_max_new);
[TFD_RGK, Phi] = rgk(Signal,2);
[TFD_FADTFD orient]=FADTFD_EEGS(Signal);
%% Plotting
figure; SetFigDef(16,9,'Times',20); tfsapl(Signal,TFD_True(:,1:tr:end),'SampleFreq',fs,'Title','(a)', 'TFfontSize' , 20,'plotfn','pcolor');
% figure; tfsapl(sig,TFD_WVD(:,1:tr:end),'SampleFreq',fs,'Title','(a)','TimePlot', 'on', 'FreqPlot', 'on', 'TFfontSize' , 20);
figure; SetFigDef(16,9,'Times',20);tfsapl(Signal,TFD_SPEC(:,1:tr:end),'SampleFreq',fs,'Title','(b)','TimePlot', 'on', 'FreqPlot', 'on', 'TFfontSize' , 20,'plotfn','pcolor');
figure; SetFigDef(16,9,'Times',20);tfsapl(Signal,TFD_EMBD(:,1:tr:end),'SampleFreq',fs,'Title','(c)','TimePlot', 'on', 'FreqPlot', 'on', 'TFfontSize' , 20,'plotfn','pcolor');
figure; SetFigDef(16,9,'Times',20);tfsapl(Signal,TFD_SM(:,1:tr:end),'SampleFreq',fs,'Title','(d)','TimePlot', 'on', 'FreqPlot', 'on', 'TFfontSize' , 20,'plotfn','pcolor');
figure; SetFigDef(16,9,'Times',20);tfsapl(Signal,TFD_CKD(:,1:tr:end),'SampleFreq',fs,'Title','(e)','TimePlot', 'on', 'FreqPlot', 'on', 'TFfontSize' , 20,'plotfn','pcolor');
figure; SetFigDef(16,9,'Times',20);tfsapl(Signal,TFD_RGK(:,1:tr:end),'SampleFreq',fs,'Title','(f)', 'TFfontSize' , 20,'plotfn','pcolor');
figure; SetFigDef(16,9,'Times',20);tfsapl(Signal,TFD_AOK(:,1:tr:end),'SampleFreq',fs,'Title','(g)','TimePlot', 'on', 'FreqPlot', 'on', 'TFfontSize' , 20,'plotfn','pcolor');
figure; SetFigDef(16,9,'Times',20);tfsapl(Signal,TFD_FADTFD(:,1:tr:end),'SampleFreq',fs,'Title','(h)','TimePlot', 'on', 'FreqPlot', 'on', 'TFfontSize' , 20,'plotfn','pcolor');
% t=0:1/32:8-1/32;
% figure;plot(t,sig);
% save simsig Signal;