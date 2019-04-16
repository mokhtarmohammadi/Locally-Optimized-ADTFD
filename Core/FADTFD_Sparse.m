function [I Id]= FADTFD_Sparse(x)
N=length(x);
alpha=0.5;
t=-floor(N/2):floor(N/2);
h1=1./((cosh(t).^2).^alpha);
h1=h1/sum(h1);
t=-floor(N/2):floor(N/2);
h2=1./((cosh(t).^2).^alpha);
h2=h2/sum(h2);
B2=h1'*h2;
S_TF= HTFD_new1(x,2,30,3*(N/4)-4);

for ii=N/4:4:N/2
    [S_TF1,orient]=HTFD_new1(x,2,30,ii);
    if (sum(sum(abs(S_TF1).^0.5)))^2/(sum(sum(S_TF1))<sum(sum(abs(S_TF).^0.5)))^2/sum(sum(S_TF)) 
        S_TF=S_TF1;
     else
         break;
    end
end

[adtfd1,orient]=HTFD_new1(x,2,30,ii-2);
adtfd1(adtfd1<0)=0;
orient=orient*3;
adtfd1(and(75<orient,orient<115))=0;
adtfd=(adtfd1);
% I=filter2(B2,adtfd);
I=adtfd;
Id=orient/3;
end