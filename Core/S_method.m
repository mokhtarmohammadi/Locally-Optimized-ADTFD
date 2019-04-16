function [ S_TF ] = S_method(tmp_sig, L,wl )
%tfd = quadtfd( tmp_sig, length(tmp_sig)-1, 1, 'specx', wL, 'hamm',256);

%[tfd, rstft, hat] = tfrrsp_sq(real(tmp_sig), 1:length(s), 2*length(s),  hamming(85), 0, 'sq');

%tfd=stft((tmp_sig),1,2,wl,1,2*length(tmp_sig));
%[tfd, rtfr, hat] = tfrrsp_sq(real(tmp_sig), 1:length(tmp_sig), 2*length(tmp_sig),  hamming(wl), 1, 'stft');
[tfd,t,f] = tfrstft(real(tmp_sig)',1:length(tmp_sig), 2*length(tmp_sig),  hamming(wl), 1);
%tfd=abs(tfd);
tfd=tfd(1:128,1:128);
P=hamming(2*L+1);
% L=1;
 P=ones(1,2*L+1);

[NN,MM]=size(tfd);
S_TF=zeros(size(tfd));

%for n=1:NN
%   for m=1:MM
for n=1:MM;%:128
    for m=1:NN
        
        for l=-L:L
            if and(and(m+l<=NN,m+l>0),and(m-l<=NN,m-l>0))
               % S_TF(n,m)=S_TF(n,m)+P(L+1+l)*tfd(n,m+l)*conj(tfd(n,m-l));
                S_TF(m,n)=S_TF(m,n)+P(L+1+l)*(tfd(m+l,n))*conj((tfd(m-l,n)));
            end
        end
    end
end
S_TF=real(S_TF);


