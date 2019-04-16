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
x=exp(1i*0.001*2*pi*n.^3/M+1i*2*2*pi*n/M)+exp(-1i*0.001*2*pi*n.^3/M+1i*45*2*pi*n/M);
%%%%%%%%%%%%%%%
        
%%%%%%%%%%%%%%%%%5
N=length(x);
iii=rand(1,M);
iii(iii<0.5)=0;
iii(iii>0)=1;
 x(iii==0)=0;
 load sig6;
 %x(1:2:end)=0;
 %x(2:3:end)=0;
%x(3:4:end)=0;
%x(3:5:end)=0;
tr=1;
% x=analytic_x;
alpha=0.5;
t=-floor(N/2):floor(N/2);
h1=1./((cosh(t).^2).^alpha);
h1=h1/sum(h1);
t=-floor(N/2):floor(N/2);
h2=1./((cosh(t).^2).^alpha);
h2=h2/sum(h2);
B2=h1'*h2;

nn=1;
for ii=0:M-1
    A(ii+1,:)=exp(1i*ii*2*pi*n/M)/M;
end

%%%%%%%%%%%
% wvd=(wvd1(x,length(x)));
fadtfd=FADTFD_Sparse(x);
[hadtfd]=HTFD_new12(x);
% S_TF  = S_method(x, 53,133 );
S_TF  = Adaptive_S_method(x, 53,'hann' );
S_TF=imresize(S_TF,[length(x) length(x)]);
I_max_new=adaptive_optimal_tfd(x);
aoktfd=real(I_max_new);
aoktfd=imresize(aoktfd,[length(x) length(x)]);
%%%%%%%%%%%%%%%%
SetFigDef(16,9,'Times',20); tfsapl(x,fadtfd(:,1:tr:128),'Title','(a)', 'TFfontSize' , 20,'plotfn','pcolor');
TF_f1=(sum(sum(abs( fadtfd).^0.5)))^2
% print('Fig6-1','-dpng','-r600');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_TF(S_TF<0)=0;
figure; SetFigDef(16,9,'Times',20); tfsapl(x,S_TF(:,1:tr:128),'Title','(b)', 'TFfontSize' , 20,'plotfn','pcolor');
TF_f2=(sum(sum(abs(S_TF).^0.5)))^2
% print('Fig6-2','-dpng','-r600');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aoktfd(aoktfd<0)=0;
figure; SetFigDef(16,9,'Times',20); tfsapl(x,aoktfd(:,1:tr:128),'Title','(c)', 'TFfontSize' , 20,'plotfn','pcolor');
TF_f3=(sum(sum(abs( aoktfd).^0.5)))^2
% print('Fig6-3','-dpng','-r600');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hadtfd(hadtfd<0)=0;
hadtfd=filter2(B2,hadtfd);
figure;SetFigDef(16,9,'Times',20); tfsapl(x,hadtfd(:,1:tr:128),'Title','(d)', 'TFfontSize' , 20,'plotfn','pcolor');
TF_f4=(sum(sum(abs( hadtfd).^0.5)))^2
% print('Fig6-4','-dpng','-r600');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cw=real(quadtfd(x,length(x)/2-1,1,'cw',64,128));
cw(cw<0)=0;
figure; SetFigDef(16,9,'Times',20); tfsapl(x,cw(:,1:tr:128),'Title','(e)', 'TFfontSize' , 20,'plotfn','pcolor');
TF_f5=(sum(sum(abs( cw).^0.5)))^2
% print('Fig6-5','-dpng','-r600');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 AFS=tfr_stft_high(hilbert(x));
 AFS(AFS<0)=0;
 figure; SetFigDef(16,9,'Times',20); tfsapl(x,AFS(:,1:tr:128),'Title','(f)', 'TFfontSize' , 20,'plotfn','pcolor');
 TF_f6=(sum(sum(abs( cw).^0.5)))^2
%  print('Fig6-6','-dpng','-r600');
%  save sig.mat x;
