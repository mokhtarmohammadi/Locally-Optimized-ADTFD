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
 load sig4;
 %x(1:2:end)=0;
 %x(2:3:end)=0;
%x(3:4:end)=0;
%x(3:5:end)=0;

% x=analytic_x;
%  x=hilbert(x);

nn=1;
for ii=0:M-1
    A(ii+1,:)=exp(1i*ii*2*pi*n/M)/M;
end
%%%%%%%%%%%%%%%%
N=length(x);
alpha=0.5;
t=-floor(N/2):floor(N/2);
h1=1./((cosh(t).^2).^alpha);
h1=h1/sum(h1);
t=-floor(N/2):floor(N/2);
h2=1./((cosh(t).^2).^alpha);
h2=h2/sum(h2);
B2=h1'*h2;

%%%%%%%%%%%
% wvd=(wvd1(x,length(x)));
FADTFD=FADTFD_Sparse(x);
[hadtfd]=HTFD_new12(x);

%%%%%%%%
%  S_TF  = S_method(x, 33,133 );
 S_TF  = Adaptive_S_method(x,53,'bart' );
 S_TF=imresize(S_TF,[length(x) length(x)]);
% g=extnd_mbd(0.5,0.5,1,length(x));
% gsig=ifft(fft(wvd.').');
% smg=gsig.*g;
% EMBD=real(fft(ifft(smg.').'));
% EMBD=imresize(EMBD,[length(x) length(x)]);
%%%%%%%
I_max_new=adaptive_optimal_tfd(x);
aoktfd=real(I_max_new);
aoktfd=imresize(aoktfd,[length(x) length(x)]);
%%%%%%%%%%%%%%%%
SetFigDef(16,9,'Times',20); tfsapl(x,FADTFD(:,1:tr:128),'Title','(a)', 'TFfontSize' , 20,'plotfn','pcolor');
TF_FADTFD=sum(sum(abs( FADTFD).^0.5))/sum(sum(FADTFD))

% print('Fig4-1','-dpng','-r600');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_TF(S_TF<0)=0;
figure; SetFigDef(16,9,'Times',20); tfsapl(x,S_TF(:,1:tr:128),'Title','(b)', 'TFfontSize' , 20,'plotfn','pcolor');
% print('Fig4-2','-dpng','-r600');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aoktfd(aoktfd<0)=0;
figure; SetFigDef(16,9,'Times',20); tfsapl(x,aoktfd(:,1:tr:128),'Title','(c)', 'TFfontSize' , 20,'plotfn','pcolor');
% print('Fig4-3','-dpng','-r600');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hadtfd(hadtfd<0)=0;
hadtfd=filter2(B2,hadtfd);
figure;SetFigDef(16,9,'Times',20); tfsapl(x,hadtfd(:,1:tr:128),'Title','(d)', 'TFfontSize' , 20,'plotfn','pcolor');
% print('Fig4-4','-dpng','-r600');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cw=real(quadtfd(x,length(x)/2-1,1,'cw',64,128));
cw(cw<0)=0;
figure; SetFigDef(16,9,'Times',20); tfsapl(x,cw(:,1:tr:128),'Title','(e)', 'TFfontSize' , 20,'plotfn','pcolor');
% print('Fig4-5','-dpng','-r600');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 AFS=tfr_stft_high(hilbert(x));
 AFS(AFS<0)=0;
 figure; SetFigDef(16,9,'Times',20); tfsapl(x,AFS(:,1:tr:128),'Title','(f)', 'TFfontSize' , 20,'plotfn','pcolor');
%  print('Fig4-6','-dpng','-r600');
%  save sig.mat x;
