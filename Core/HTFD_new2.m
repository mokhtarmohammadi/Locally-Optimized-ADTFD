% Please cite "Mokhtar Mohammadi, Ali Akbar Pouyan, Nabeel Ali Khan, Vahid Abolghasemi:
%              Locally Optimized Adaptive Directional Time-Frequency Distributions. 
%               CSSP 37(8): 3154-3174 (2018)" 

function [Inew2 Id2]= HTFD_new2(s,alpha1,alpha2,N)

%% parameter


M=180;
% R=1;
R=3;
times=1;
alpha=0.3;
addpath('D:\classes\time-frequency\tfsa_6-0\tfsa_6-0\windows\win64_bin'); % utilisation de la toolbox tfsa

%% EMBD
I = quadtfd(s, length(s)-1, 1, 'wvd',length(s)*times);

t=-floor(N/2):floor(N/2);
h1=1./((cosh(t).^2).^alpha);
h1=h1/sum(h1);
t=-floor(N/2):floor(N/2);
h2=1./((cosh(t).^2).^alpha);
h2=h2/sum(h2);
B2=h1'*h2;
 I=filter2(B2,I); % filtrage EMBD
% I=tfr_stft_high(hilbert(s));
% I=filter2(B2,I);
%% Mask construction
[X,Y]=meshgrid(-1:2/N:1,-1:2/N:1); % discretisation de l'espace
%[X,Y]=meshgrid(-2:4/N:2,-2:4/N:2); % discretisation de l'espace

%% Init

PQ = paddedsize(size(I));
F = fft2(double(I), PQ(1), PQ(2));
FIabs = fft2(double(abs(I)), PQ(1), PQ(2));


A=exp((-1/2)*(((alpha1*X).^2)+(alpha2*Y).^2)); %gaussian mask creation
A=A.*(1-alpha2*alpha2*Y.^2);
A=A/sum(sum(abs(A))); %normalization
%mesh(A)
H = fft2(double(A), PQ(1), PQ(2));

F_fH = H.*FIabs;
ffi = abs(ifft2(F_fH));
IImax = (ffi(round(length(A)/2):end-round(length(A)/2), round(length(A)/2):end-round(length(A)/2)));
positionmax = ones(size(IImax));

F_fH = H.*F;
ffi = (real(ifft2(F_fH)));
Inew = ffi(round(length(A)/2):end-round(length(A)/2), round(length(A)/2):end-round(length(A)/2));

%% filtering in the ambiguity domain
for ii=0:M/R-1
    
    angle=pi*ii*R/M; % degres->radian
    
    X1=X*cos(angle)-Y*sin(angle); %axis rotation
    Y1=X*sin(angle)+Y*cos(angle); %axis rotation
    A=exp((-1/2)*(((alpha1*X1).^2)+(alpha2*Y1).^2)); %gaussian mask creation
    A=A.*(1-alpha2*alpha2*Y1.^2);
    A=A/sum(sum(abs(A))); %normalization
    
    H = fft2(double(A), PQ(1), PQ(2));
    F_fH = H.*FIabs;
    ffi = (ifft2(F_fH));
    ffi=abs(ffi);
    IItemp = (ffi(round(length(A)/2):end-round(length(A)/2), round(length(A)/2):end-round(length(A)/2)));
    IImax_old = IImax;
    IImax=max(IImax,IItemp);
%     IImaxc=(repcolumn((mean(IImax))',MM))';
%     IImaxr=repcolumn(mean(IImax,2),MM);

    positionChangement = (IImax-IImax_old)~=0;% & IImax> IImaxc & IImax> IImaxr;
    positionmax(positionChangement) = (ii+1)*3;
    
    F_fH = H.*F; %
    ffi = (real(ifft2(F_fH)));
    II3temp = ffi(round(length(A)/2):end-round(length(A)/2), round(length(A)/2):end-round(length(A)/2));
    Inew(positionChangement)= II3temp(positionChangement);
end

[M1,N1]=size(I);
Id=positionmax;
Id2=Id(1:M1,1:N1);
Inew2 = Inew(1:M1,1:N1);
Inew2=filter2(B2,Inew2);
Inew2(Inew2<0)=0;
%Inew2=I;
%plotTF(Inew2,1);
end