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

M=128;
n=0:M-1;
 x=exp(1i*0.0007*2*pi*n.^3/M+1i*2*2*pi*n/M)+exp(1i*0.0007*2*pi*n.^3/M+1i*16*2*pi*n/M);
%%%%%%%%%%%%%%%
        tr=1;
%%%%%%%%%%%%%%%%%5
N=length(x);
iii=rand(1,M);
iii(iii<0.5)=0;
iii(iii>0)=1;
  x(iii==0)=0;
 load sig6;
 
nn=1;
for ii=0:M-1
    A(ii+1,:)=exp(1i*ii*2*pi*n/M)/M;
end
%%%%%%%%%%%%%%%%
N=length(x);
%%%%%%%%%%%
% wvd=(wvd1(x,length(x)));
FADTFD=FADTFD_Sparse(x);
[hadtfd]=HTFD_new12(x);

%%%%%%%%
 S_TF  = Adaptive_S_method(x,53,'hann' );
 S_TF=imresize(S_TF,[length(x) length(x)]);
 I_max_new=adaptive_optimal_tfd(x);
tfdAOK=real(I_max_new);
tfdAOK=imresize(tfdAOK,[length(x) length(x)]);
%%%%%%%%%%%%%%%%
TF_FADTFD=(sum(sum(abs( FADTFD).^0.5)))^2/sum(sum(FADTFD));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_TF(S_TF<0)=0;
TF_A_S=(sum(sum(abs( S_TF).^0.5)))^2/sum(sum(S_TF));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tfdAOK(tfdAOK<0)=0;
TF_AOK=(sum(sum(abs( tfdAOK).^0.5)))^2/sum(sum(tfdAOK));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hadtfd(hadtfd<0)=0;
TF_HADTFD=(sum(sum(abs( hadtfd).^0.5)))^2/sum(sum(hadtfd));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cw=real(quadtfd(x,length(x)/2-1,1,'cw',64,128));
cw(cw<0)=0;
TF_CW=(sum(sum(abs( cw).^0.5)))^2/sum(sum(cw));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 AFS=tfr_stft_high(hilbert(x));
 AFS(AFS<0)=0;
TF_AFS=(sum(sum(abs( AFS).^0.5)))^2/sum(sum(AFS));
%  save sig.mat x;
fprintf('AFS  = %0.f\n', TF_AFS/100);
fprintf('AOK  = %0.f\n', TF_AOK/100);
fprintf('HADTFD  = %0.f\n', TF_HADTFD/100);
fprintf('FADTFD  = %0.f\n', TF_FADTFD/100);
fprintf('Adaptive S_method  = %0.f\n', TF_A_S/100);
fprintf('Choi-William  = %0.f\n', TF_CW/100);
