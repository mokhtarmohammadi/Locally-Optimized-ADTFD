% Please cite "Mokhtar Mohammadi, Ali Akbar Pouyan, Nabeel Ali Khan, Vahid Abolghasemi:
%              Locally Optimized Adaptive Directional Time-Frequency Distributions. 
%               CSSP 37(8): 3154-3174 (2018)" 

function [Inew2 k Id] = HTFD_AD2(s)

M=180;
% R=1;
R=3;
times=1;

%%%%%%%%%%%%
%% filtering in the ambiguity domain
% degres->radian

i=1;
%%%%%%%%Optimized WL based on the measure
L=[   32    64                            32       64    ];
alpha1= [     3     3                           2        2            ];
alpha2= [     6      8                         20         30        ];  

for i=1:length(L)
    
[Inew2(i,:,:) Id1(i,:,:)]=HTFD_new1(s,alpha1(i),alpha2(i),L(i));
end
k=length(L);
Id=Id1;