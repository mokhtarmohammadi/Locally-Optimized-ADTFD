%  
%%  Author:
%     Mokhtar Mohammadi
%In this code we assume that the user add in this path the TFSAP toolboox.
% Please cite "Mokhtar Mohammadi, Ali Akbar Pouyan, Nabeel Ali Khan, Vahid Abolghasemi:
%              Locally Optimized Adaptive Directional Time-Frequency Distributions. 
%               CSSP 37(8): 3154-3174 (2018)" 

%% Locating and Adding Functions Directory
clear all;
currentDirectory = pwd;
[upperPath, ~, ~] = fileparts(currentDirectory);
addpath([upperPath '\Core'])
addpath([upperPath '\TFSA5'])
addpath([upperPath '\Core\AOK'])
addpath([upperPath '\TFTB2'])

t=0:255;
s1=cos(0.1*pi.*t +0.000002*pi.*(t.^3));
s2=cos(0.05*pi.*t +0.000002*pi.*(t.^3));
s3=cos(0.9*pi.*t -0.000001*pi.*(t.^3));
s4=cos(0.95*pi.*t -0.000001*pi.*(t.^3));
s5=2*exp(-0.0015*(t-55).^2).*cos(0.6.*t*pi);
s6=2*exp(-0.0015*(t-180).^2).*cos(0.6.*t*pi);
x=(s1+s2+0.5*s3+0.5*s4+s5+s6);
% s=1*cos(2*pi*(0.025*t+0.000001*t.^3))+1*cos(2*pi*(0.05*t+0.000001*t.^3))+0.5*cos(2*pi*(0.45*t-0.0000005*t.^3))+0.5*cos(2*pi*(0.475*t-0.0000005*t.^3));
% x=s;
% N_C=4;
display('This code shows TFDs in the following order');
display('1. AFS');
display('2. AOK-TFD');
display('3. Adaptive directional TFD (3,8)');
display('4. HADTFD');
display('5. FATFD');
display('6. Spectrogram');
display('7. SPWVD');
display('8. Reasigned SPWVD');
%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%% AFS %%%%%%%%
SetFigDef(16,9,'Times',20);
AFS=tfr_stft_high(hilbert(x));
tfsapl(x,AFS, 'TFfontSize' , 20,'Title','(a)','plotfn','pcolor');
title('Adaptive fractional Spectrogram');
TF_f1m=sum(sum(abs( AFS).^0.5)) /sum(sum(AFS))%% The measure for Table1
TF_f1s=(sum(sum(abs( AFS).^0.5)))^2/sum(sum(AFS))
% print('Fig2-1','-dpng','-r600');
%%%%%%%%%%%%% AOK TFD%%%%%%%%%
I_max_new=adaptive_optimal_tfd(x);
tfdAOK=real(I_max_new);
figure;SetFigDef(16,9,'Times',20);
tfsapl(x,tfdAOK, 'TFfontSize' , 20,'Title','(b)','plotfn','pcolor');
title('AOK-TFD');
TF_f2m=sum(sum(abs( tfdAOK).^0.5)) /sum(sum(tfdAOK))
TF_f2s=(sum(sum(abs( tfdAOK).^0.5)))^2/sum(sum(tfdAOK))
% print('Fig2-2','-dpng','-r600');
% %%%%%%%%%%%%%%%% ADTFD(3,8,32) %%%%%%%%%%%
ADTFD2 = HTFD_new2(x, 3,8,32);
ADTFD2(ADTFD2<0)=0;
ADTFD2=imresize(ADTFD2,[length(x) length(x)]);
figure; SetFigDef(16,9,'Times',20); tfsapl(x,ADTFD2, 'TFfontSize' , 20,'Title','(c)','plotfn','pcolor');
title('Adaptive direction TFD');
TF_f3m=sum(sum(abs( ADTFD2).^0.5)) /sum(sum(ADTFD2))%% The measure for Table1
TF_f3s=(sum(sum(abs( ADTFD2).^0.5)))^2
% print('Fig2-3','-dpng','-r600');
%%%%%%%%%%%%%%%%%% HADTFD %%%%
[HADTFD] = HTFD_new12(x);
HADTFD(HADTFD<0)=0;
HADTFD=imresize(HADTFD,[length(x) length(x)]);
figure; SetFigDef(16,9,'Times',20); tfsapl(x,HADTFD, 'TFfontSize' , 20,'Title','(d)','plotfn','pcolor');
title(' HADTFD');
TF_f4m=sum(sum(abs( HADTFD).^0.5)) /sum(sum(HADTFD))%% The measure for Table1
TF_f4s=(sum(sum(abs( HADTFD).^0.5)))^2
% print('Fig2-4','-dpng','-r600');
%%%%%%%%%%%%%%%%%% FADTFD %%%
FADTFD= HTFD_new10(x);
FADTFD(FADTFD<0)=0;
FADTFD=imresize(FADTFD,[length(x) length(x)]);
figure; SetFigDef(16,9,'Times',20); tfsapl(x,FADTFD, 'TFfontSize' , 20,'Title','(e)','plotfn','pcolor');
title('Fully automatic adaptive TFD new idea');
TF_f5m=sum(sum(abs( FADTFD).^0.5)) /sum(sum(FADTFD))%% The measure for Table1
TF_f5=(sum(sum(abs( FADTFD).^0.5)))^2/sum(sum(FADTFD))
% print('Fig2-5','-dpng','-r600');
%%%%%%%%%%%%%%%%%%% Spectrogram %%%%%%%%%
spec = quadtfd( x, length(x)-1, 1, 'specx', 85, 'hamm',256);
figure; SetFigDef(16,9,'Times',20); tfsapl(x,spec, 'TFfontSize' , 20,'Title','(f)','plotfn','pcolor');
title(' Spectrogram');
TF_f6m=sum(sum(abs( spec).^0.5)) /sum(sum(spec))%% The measure for Table1
TF_f6s=(sum(sum(abs( spec).^0.5)))^2
% print('Fig2-6','-dpng','-r600');
%%%%%%%%%%%%%%%%%SPWVD and its Reasigend version%%%%%
[tfr,rtfr]=tfrrspwv(hilbert(x'));
tfr(tfr<0)=0;
rtfr(rtfr<0)=0;
figure; SetFigDef(16,9,'Times',20); tfsapl(x',tfr, 'TFfontSize' , 20,'Title','(g)','plotfn','pcolor');
% print('Fig2-7','-dpng','-r600');
figure; SetFigDef(16,9,'Times',20); tfsapl(x',rtfr, 'TFfontSize' , 20,'Title','(h)','plotfn','pcolor');
%print('Fig2-8','-dpng','-r600');
TF_f7m=sum(sum(abs( tfr).^0.5)) /sum(sum(tfr))%% The measure for Table1
TF_f7s=(sum(sum(abs( tfr).^0.5)))^2
TF_f8m=sum(sum(abs( rtfr).^0.5)) /sum(sum(rtfr))%% The measure for Table1
TF_f8s=(sum(sum(abs( rtfr).^0.5)))^2

