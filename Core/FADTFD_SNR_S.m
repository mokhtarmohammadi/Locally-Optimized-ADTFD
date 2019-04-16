function I_max_new = FADTFD_SNR_S(x)
N=length(x);
alpha=0.5;
t=-floor(N/2):floor(N/2);
h1=1./((cosh(t).^2).^alpha);
h1=h1/sum(h1);
t=-floor(N/2):floor(N/2);
h2=1./((cosh(t).^2).^alpha);
h2=h2/sum(h2);
B2=h1'*h2;

[fadtfd,orient]= HTFD_new1(x, 2,30,32);
fadtfd(fadtfd<0)=0;
orient=orient*3;
fadtfd(and(75<orient,orient<115))=0;
fadtfd=filter2(B2,fadtfd);
I_max_new=fadtfd;
I_max_new=imresize(I_max_new,[length(x) length(x)]);
