% Please cite "Mokhtar Mohammadi, Ali Akbar Pouyan, Nabeel Ali Khan, Vahid Abolghasemi:
%              Locally Optimized Adaptive Directional Time-Frequency Distributions. 
%               CSSP 37(8): 3154-3174 (2018)" 

function [I Id]= HTFD_new10(s)

MM=32;
M=180;
% R=1;
R=3;
times=1;

%addpath('D:\classes\time-frequency\tfsa_6-0\tfsa_6-0\windows\win64_bin'); % utilisation de la toolbox tfsa
N=128;
alpha=0.5;
I = quadtfd(s, length(s)-1, 1, 'wvd',length(s)*times);

t=-floor(N/2):floor(N/2);
h1=1./((cosh(t).^2).^alpha);
h1=h1/sum(h1);
t=-floor(N/2):floor(N/2);
h2=1./((cosh(t).^2).^alpha);
h2=h2/sum(h2);
B2=h1'*h2;


I = quadtfd(s, length(s)-1, 1, 'wvd',length(s)*times);
[MM,NN]=size(I);
%%%%%%%%%%%%%%%%%%%%%%%%%
[ALLabstfr k ID]=HTFD_AD2all(s);
II=ALLabstfr(1,:,:);
Id=ID(1,:,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
position=zeros(size(I));
for    i=2:k
     %   Itemp=  ALLabstfr(i,:,:);

   % figure; tfsapl(s,squeeze(II));
%     positionChangement =  (sum(sum(abs( squeeze(ALLabstfr(i,:,:))).^0.5)))^2<(sum(sum(abs( squeeze(II)).^0.5)))^2;
tfd=squeeze(ALLabstfr(i,:,:));
tfd1=squeeze(II);
% tfdc=(repcolumn((mean(tfd))',MM))';
% tfdr=repcolumn(mean(tfd,2),MM);
[x,y]=size(tfd);
I1=zeros(x+MM,y+MM);
I1(MM/2+1:end-MM/2,MM/2+1:end-MM/2)=tfd;
for ii=MM/2+1:length(s)+MM/2
    for jj=MM/2+1:length(s)+MM/2

        B=I1(ii-MM/2:ii+MM/2,jj-MM/2:jj+MM/2);
        tfdc(ii-MM/2,jj-MM/2)=sum(I1(ii-MM/2:ii+MM/2,jj));
        tfdr(ii-MM/2,jj-MM/2)=sum(I1(ii,jj-MM/2:jj+MM/2));
    end
end

% tfd1c=(repcolumn((mean(tfd1))',MM))';
% tfd1r=repcolumn(mean(tfd1,2),MM);

positionChangement = tfd<tfd1;%(tfd<tfd1| tfd>tfdc & tfd>tfdr);% & (w.NumObjects < w1.NumObjects);
II(positionChangement)=ALLabstfr(i,positionChangement);
Id(positionChangement)= ID(i, positionChangement);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[M1,N1]=size(I);
I = squeeze(II(1,1:M1,1:N1));
Id = squeeze(Id(1,1:M1,1:N1));
%I=AEMBD1(I);
I=filter2(B2,I);
I(I<0)=0;


% figure(2);
% tfsapl(s,I);
%plotTFS(I',1);

end