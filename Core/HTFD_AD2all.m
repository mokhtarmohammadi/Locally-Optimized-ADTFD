function [Inew2 k Id] = HTFD_AD2all(s)

M=180;
% R=1;
R=3;
times=1;

%%%%%%%%%%%%
%% filtering in the ambiguity domain
% degres->radian

i=1;
%%%%%%%%Optimized WL based on the measure
L=[  32 32    64         32     54        54  32    32 64    64   64    64 64];
alpha1= [ 2    3     3     3     3          3    2    2     2 2   2    2   2     ];
alpha2= [   2  6      8      4     7         5    19  20    23  27   30  25 15     ];  

for i=1:length(L)
    
[Inew2(i,:,:) Id1(i,:,:)]=HTFD_new1(s,alpha1(i),alpha2(i),L(i));
end
k=length(L);
Id=Id1;