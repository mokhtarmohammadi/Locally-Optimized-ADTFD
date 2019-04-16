function [ S_TF ] = Adaptiv_S_method(Signal,wl,win)

TFD_SPEC = quadtfd( Signal, length(Signal)-1, 1, 'specx', wl, win,256);
N=size(TFD_SPEC,1);
 S_TF=SM_calc(TFD_SPEC,0);
for L=1:N
    S_TF1=SM_calc(TFD_SPEC,L);
       if (sum(sum(abs(S_TF1).^0.5)))^2/sum(sum(S_TF1))<(sum(sum(abs(S_TF).^0.5)))^2/sum(sum(S_TF))
      S_TF=S_TF1;
       else
           break;
end
end
