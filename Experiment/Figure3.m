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

t=0:255;
s1=cos(0.1*pi.*t +0.000003*pi.*(t.^3));
s2=cos(0.65*pi.*t -0.000003*pi.*(t.^3));
s3=2*exp(-0.0015*(t-94).^2).*cos(0.8.*t*pi);
s4=2*exp(-0.0015*(t-94).^2).*cos(0.4.*t*pi);
x=(s1+s2+s3+s4);

N_C=4;
display('This code shows TFDs in the following order');
display('1. AFS');
display('2. AOK-TFD');
display('3. Adaptive directional TFD (3,8)');
display('4. HADTFD');
display('5. FATFD');
display('6. Spectrogram');
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% AFS %%%%%%%%
AFS=tfr_stft_high(hilbert(x));
SetFigDef(16,9,'Times',20); tfsapl(x,AFS, 'TFfontSize' , 20,'Title','(a)','plotfn','pcolor');
title('Adaptive fractional Spectrogram');
% print('Fig3-1','-dpng','-r600');
%%%%%%%%%%%%% AOK TFD%%%%%%%%%
I_max_new=adaptive_optimal_tfd(x);
tfdAOK=real(I_max_new);
figure; SetFigDef(16,9,'Times',20); tfsapl(x,tfdAOK, 'TFfontSize' , 20,'Title','(b)','plotfn','pcolor');
title('AOK-TFD');
% print('Fig3-2','-dpng','-r600');
% %%%%%%%%%%%%%%%% ADTFD(3,8,32) %%%%%%%%%%%
ADTFD2 = HTFD_new2(x, 2,20,102);
ADTFD2(ADTFD2<0)=0;
ADTFD2=imresize(ADTFD2,[length(x) length(x)]);
figure; tfsapl(x,ADTFD2, 'TFfontSize' , 20,'Title','(c)','plotfn','pcolor');
title('Adaptive direction TFD');
% print('Fig3-3','-dpng','-r600');
%%%%%%%%%%%%%%%%%% HADTFD %%%%
[HADTFD] = HTFD_new12(x);
HADTFD(HADTFD<0)=0;
HADTFD=imresize(HADTFD,[length(x) length(x)]);
figure; tfsapl(x,HADTFD, 'TFfontSize' , 20,'Title','(d)','plotfn','pcolor');
title(' HADTFD');
% print('Fig3-4','-dpng','-r600');
%%%%%%%%%%%%%%%%%% FADTFD %%%
FADTFD= HTFD_new10(x);
FADTFD(FADTFD<0)=0;
FADTFD=imresize(FADTFD,[length(x) length(x)]);
figure; SetFigDef(16,9,'Times',20); tfsapl(x,FADTFD, 'TFfontSize' , 20,'Title','(e)','plotfn','pcolor');
title('Fully automatic adaptive TFD new idea');
% print('Fig3-5','-dpng','-r600');
%%%%%%%%%%%%%%%%%% Spectrogram %%%%%%%%%
spec = quadtfd( x, length(x)-1, 1, 'specx', 85, 'hamm',256);
figure; SetFigDef(16,9,'Times',20); tfsapl(x,spec, 'TFfontSize' , 20,'Title','(f)','plotfn','pcolor');
title(' Spectrogram');
% print('Fig3-6','-dpng','-r600');


