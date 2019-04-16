% Please cite "Mokhtar Mohammadi, Ali Akbar Pouyan, Nabeel Ali Khan, Vahid Abolghasemi:
%              Locally Optimized Adaptive Directional Time-Frequency Distributions. 
%               CSSP 37(8): 3154-3174 (2018)" 

function [TF_EEG orient]= FADTFD_EEGR(x)

%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(x);
S_TF= HTFD_new1(x,2,30,3*(N/8));
i_opt=3*(N/8);
for ii=3*(N/8):4:N/2
    [S_TF1,orient]=HTFD_new1(x,2,30,ii);
    if (sum(sum(abs(S_TF1).^0.5)))^2/sum(sum(S_TF1))<(sum(sum(abs(S_TF).^0.5)))^2/sum(sum(S_TF)) 
        S_TF=S_TF1;
     else
         break;
    end
end
    
[TFD_ADTFD1 orient1]=F_ADTFD(x,2,30,ii-2,[5 175]);
Orient1=orient1*3;
  [TFD_ADTFD2 orient2]=F_ADTFD(x,3,8,42,[95 85]);
Orient2=orient2*3;
        
 TFD_ADTFD=(TFD_ADTFD1+TFD_ADTFD2)/2;
 Orient=(Orient1+Orient2)/2;
 TFD_ADTFD(and(75<Orient,Orient<115))=0;
 orient=orient1;
 TF_EEG=TFD_ADTFD;
end
