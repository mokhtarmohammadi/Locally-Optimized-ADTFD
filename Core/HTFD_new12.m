% Please cite "Mokhtar Mohammadi, Ali Akbar Pouyan, Nabeel Ali Khan, Vahid Abolghasemi:
%              Locally Optimized Adaptive Directional Time-Frequency Distributions. 
%               CSSP 37(8): 3154-3174 (2018)" 

function I= HTFD_new12(s)

%% parameter
N=256;
M=180;
% R=1;
R=3;
times=1;
addpath('D:\classes\time-frequency\tfsa_6-0\tfsa_6-0\windows\win64_bin'); % utilisation de la toolbox tfsa
 I = quadtfd(s, length(s)-1, 1, 'wvd',length(s)*times);
%%%%%%%%%%%%%%%%%%%%%%%%%
[M1 N1]=size(I);
[ALLabstfr k]=HTFD_AD3(s);
%%%%%%%%%%%%%%%%%%%%%%%%
diffMem = zeros(k-2,M1,N1);
for i=1:k-2
    diff = squeeze(ALLabstfr(i,:,:)) - squeeze(ALLabstfr(i+1,:,:));
    diffMem(i,:,:) = diff;
end
clear diff;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Method of point selection
[~,Imin] = min(abs(diffMem));
clear diffMem; 

NEWabstfr = zeros(size(Imin,2),size(Imin,3));
for i=1:size(Imin,2)
    for j=1:size(Imin,3)
        NEWabstfr(i,j) = ALLabstfr(Imin(1,i,j),i,j);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [M1,N1]=size(I);
 I = NEWabstfr(1:M1,1:N1);
 I(I<0)=0;

end