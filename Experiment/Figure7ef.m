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
addpath([upperPath '\TFSA5'])
addpath([upperPath '\Core\AOK'])
addpath([upperPath '\TFTB2'])
%%%%%%%%%%%%
type=9;
[t,x,fs,IF_O]=signal_type_new(type);
x=x';
 x=awgn(x,10,'measured');
 load sig7;
% save sig7 x;
%%%%%%%%%%%% AOK TFD%%%%%%%%%
I_max_new=adaptive_optimal_tfd(x);
tfdAOK=real(I_max_new);
figure; SetFigDef(16,9,'Times',20); tfsapl(x,tfdAOK, 'TFfontSize' , 20,'Title','(b)','plotfn','pcolor');
title('AOK-TFD');
%%%%%%%%%%%%%%%%%%FADTFD%%%%%%%%%%%%%
% FADTFD= HTFD_new10(x);
FADTFD= FADTFD_SNR_C(x);
figure; SetFigDef(16,9,'Times',20); tfsapl(x,FADTFD, 'TFfontSize' , 20,'Title','(e)','plotfn','pcolor');
title('Fully automatic adaptive TFD new idea');
% save sig7 x;