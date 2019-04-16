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
s1=cos(0.1*pi.*t +0.000002*pi.*(t.^3));
s2=cos(0.05*pi.*t +0.000002*pi.*(t.^3));
s3=cos(0.9*pi.*t -0.000001*pi.*(t.^3));
s4=cos(0.95*pi.*t -0.000001*pi.*(t.^3));
s5=2*exp(-0.0015*(t-55).^2).*cos(0.6.*t*pi);
s6=2*exp(-0.0015*(t-180).^2).*cos(0.6.*t*pi);
x=(s1+s2+0.5*s3+0.5*s4+s5+s6);
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
AFS=tfr_stft_high(hilbert(x));
TF_AFS=(sum(sum(abs( AFS).^0.5)))^2/sum(sum(AFS));
% print('Fig2-1','-dpng','-r600');
%%%%%%%%%%%%% AOK TFD%%%%%%%%%
I_max_new=adaptive_optimal_tfd(x);
tfdAOK=real(I_max_new);
TF_AOK=(sum(sum(abs( tfdAOK).^0.5)))^2/sum(sum(tfdAOK));
% %%%%%%%%%%%%%%%% ADTFD(3,8,32) %%%%%%%%%%%
ADTFD2 = HTFD_new2(x, 3,8,32);
ADTFD2(ADTFD2<0)=0;
ADTFD2=imresize(ADTFD2,[length(x) length(x)]);
TF_ADTFD2=(sum(sum(abs( ADTFD2).^0.5)))^2/sum(sum(ADTFD2));
%%%%%%%%%%%%%%%%%% HADTFD %%%%
[HADTFD] = HTFD_new12(x);
HADTFD(HADTFD<0)=0;
HADTFD=imresize(HADTFD,[length(x) length(x)]);
TF_HADTFD=(sum(sum(abs( HADTFD).^0.5)))^2/sum(sum(HADTFD));
%%%%%%%%%%%%%%%%%% FADTFD %%%
FADTFD= HTFD_new10(x);
FADTFD(FADTFD<0)=0;
FADTFD=imresize(FADTFD,[length(x) length(x)]);
TF_FADTFD=(sum(sum(abs( FADTFD).^0.5)))^2/sum(sum(FADTFD));
%%%%%%%%%%%%%%%%%%% Spectrogram %%%%%%%%%
TF_SPEC = quadtfd( x, length(x)-1, 1, 'specx', 85, 'hamm',256);
TF_SPEC=(sum(sum(abs( TF_SPEC).^0.5)))^2/sum(sum(TF_SPEC));
%%%%%%%%%%%%%%%%%SPWVD and its Reasigend version%%%%%
[TF_TFR,TF_RTFR]=tfrrspwv(hilbert(x'));
TF_TFR(TF_TFR<0)=0;
TF_RTFR(TF_RTFR<0)=0;
TF_TFR=(sum(sum(abs( TF_TFR).^0.5)))^2 /sum(sum(TF_TFR));%% The measure for Table1
TF_RTFR=(sum(sum(abs( TF_RTFR).^0.5)))^2 /sum(sum(TF_RTFR));%% The measure for Table1
%%%%%%%%%%%%%%%%%%%%%%%
fprintf('AFS  = %0.f\n', TF_AFS/100);
fprintf('AOK  = %0.f\n', TF_AOK/100);
fprintf('ADTFD2  = %0.f\n', TF_ADTFD2/100);
fprintf('HADTFD  = %0.f\n', TF_HADTFD/100);
fprintf('FADTFD  = %0.f\n', TF_FADTFD/100);
fprintf('SPECTROGRAM  = %0.\n', TF_SPEC/100);
fprintf('SPWVD  = %0.f\n', TF_TFR/100);
fprintf('RSPWVD  = %0.f\n', TF_RTFR/100);
