function SM=SM_calc(S,L)
%The S-method calculation with L correction terms 
%(L=0 no correction, L=N/2-1 is the WD, all correction)
%Input S is the STFT matrix (frequency index is row index in matrix)
N=size(S,1);
SM=abs(S).^2;
for k=1:L;
    SM(1+k:N-k,:)=SM(1+k:N-k,:)+2*real(S(1:N-2*k,:).*conj(S(1+2*k:N,:)));
end
end
