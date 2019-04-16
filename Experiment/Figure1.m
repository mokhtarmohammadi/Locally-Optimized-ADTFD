%  
%%  Author:
%     Mokhtar Mohammadi
% In this code we assume that the user add in this path the TFSAP toolboox.

% Please cite "Mokhtar Mohammadi, Ali Akbar Pouyan, Nabeel Ali Khan, Vahid Abolghasemi:
%              Locally Optimized Adaptive Directional Time-Frequency Distributions. 
%               CSSP 37(8): 3154-3174 (2018)" 
%% Locating and Adding Functions Directory
clear all;
currentDirectory = pwd;
[upperPath, ~, ~] = fileparts(currentDirectory);
addpath([upperPath '\Core'])
addpath([upperPath '\TFSA5'])
t=0:255;
s1=cos(0.1*pi.*t +0.000002*pi.*(t.^3));
s2=cos(0.05*pi.*t +0.000002*pi.*(t.^3));
s3=cos(0.90*pi.*t -0.000001*pi.*(t.^3));
s4=cos(0.95*pi.*t -0.000001*pi.*(t.^3));
x=s1+s2+s3+s4;

N=length(x);
% N_C=4;
display('This code shows TFDs in the following order');
display('1. Adaptive directional TFD (2,30,16)');
display('2. Adaptive directional TFD (2,30,optimom)');
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% ADTFD(2,30,16) %%%%%%%%%%%
ADTFD1 = HTFD_new2(x, 2,30,16);
ADTFD1(ADTFD1<0)=0;
ADTFD1=imresize(ADTFD1,[length(x) length(x)]);
figure;SetFigDef(16,9,'Times',20);
tfsapl(x,ADTFD1, 'TFfontSize' , 20,'Title','(a)','plotfn','pcolor');
% print('Fig1-1','-dpng','-r600');
title('Adaptive direction TFD');
TF_f1=(sum(sum(abs( ADTFD1).^0.5)))^2
%%%%%%%%%%%%%%%%%% ADTFD with optimal WL %%%%
ADTFD1 = HTFD_new2(x, 2,30,94);
ADTFD1(ADTFD1<0)=0;
ADTFD1=imresize(ADTFD1,[length(x) length(x)]);
figure;SetFigDef(16,9,'Times',20);
tfsapl(x,ADTFD1, 'TFfontSize' , 20,'Title','(a)','plotfn','pcolor');
% print('Fig1-1','-dpng','-r600');
title('Adaptive direction TFD');
TF_f2=(sum(sum(abs( ADTFD1).^0.5)))^2
