%  
%%  Author:
%     Mokhtar Mohammadi
%In this code we assume that the user add in this path the TFSAP toolboox.
% Please cite "Mokhtar Mohammadi, Ali Akbar Pouyan, Nabeel Ali Khan, Vahid Abolghasemi:
%              Locally Optimized Adaptive Directional Time-Frequency Distributions. 
%               CSSP 37(8): 3154-3174 (2018)" 

%% Locating and Adding Functions Directory
clear all;
currentDirectory = pwd;
[upperPath, ~, ~] = fileparts(currentDirectory);
addpath([upperPath '\Core'])
addpath([upperPath '\TFSA5'])
addpath([upperPath '\Core\AOK'])
addpath([upperPath '\TFTB2'])

dis=0;
lll=1;
LL=3;
N_S=50;
type=8;
[t,x,fs,IF_O]=signal_type_new(type);
N_C=4;

xn=x;
for snr=-5:2:16
    
    for k1=1:N_S
        [t,x,fs,IF_O]=signal_type_new(type);
        
        x=awgn(x,snr,'measured',lll*200+k1+3+10);
        I_max_new=tfr_stft_high(hilbert(x'));
        IF_COMPUTE_new_edge_link;
        
    end
    var_snr_AFS(lll,:)=(mean(mse1))';
    for k1=1:N_S
        [t,x,fs,IF_O]=signal_type_new(type);
        
        x=awgn(x,snr,'measured',lll*200+k1+3+10);
        
        
        I_max_new=FADTFD_SNR_C(x);
        IF_COMPUTE_new_edge_link;
        %IF_COMPUTE;
        
    end
    var_snr_HTFD(lll,:)=(mean(mse1))';
    for k1=1:N_S
        [t,x,fs,IF_O]=signal_type_new(type);
        
        x=awgn(x,snr,'measured',lll*200+k1+3+10);
        
        adaptive_optimal_tfd;
        I_max_new=real(I_max_new);
        I_max_new(I_max_new<0)=0;
        IF_COMPUTE_new_edge_link
        %IF_COMPUTE;
        
    end
    var_snr_AOK(lll,:)=(mean(mse1))';
    lll=lll+1;
end


snr=-5:2:16
for i=1:4
    figure;
    plot(snr,10*(log10(var_snr_AOK(:,i))),'-co','linewidth',4);
    
    
    hold on;
    plot(snr, 10*(log10(var_snr_AFS(:,i))),'-rh','linewidth',4);
    
    hold on;
    plot(snr, 10*(log10(var_snr_HTFD(:,i))),'-.k+','linewidth',4);
    
    
    legend(' AOK','AFS',' FADTFD','Location','Best')
    
    xlabel('Signal to Noise Ratio');
    ylabel('Mean Square Error (dB)');
end