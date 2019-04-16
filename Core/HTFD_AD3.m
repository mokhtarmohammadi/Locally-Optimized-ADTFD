function [Inew2 k Id] = HTFD_AD3(s)

M=180;
% R=1;
R=3;
times=1;

%%%%%%%%%%%%
%% filtering in the ambiguity domain
% degres->radian

i=1;
%  L=[   12    14                            20       38    ];
% alpha1= [     3     13                           12        2            ];
% alpha2= [     16      8                         6         11        ];  
% %   L=[  30  34    60                            84      94  120    ];
% % alpha1= [    3   3     3                           2                2    2];
% alpha2= [    4   6      8                         20               30  35]; 
L=[   98    88                            22       12   ];
alpha1= [     3     3                           2        2            ];
alpha2= [     6      8                         20         30         ];  

for i=1:length(L)
    
[Inew2(i,:,:) Id1(i,:,:)]=HTFD_new1(s,alpha1(i),alpha2(i),L(i));
end
k=length(L);
Id=Id1;