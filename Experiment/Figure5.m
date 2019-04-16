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
type=2;
%type=9;
[t,x,fs,IF_O]=signal_type_new(type);
N_C=2;
xn=x;
N=length(x);
alpha=0.5;
t=-floor(N/2):floor(N/2);
h1=1./((cosh(t).^2).^alpha);
h1=h1/sum(h1);
t=-floor(N/2):floor(N/2);
h2=1./((cosh(t).^2).^alpha);
h2=h2/sum(h2);
B2=h1'*h2;
%type=15
%for snr=-5:2:16
%    for snr=-5:2:16
 
 for snr=-5:2:16
     for k1=1:N_S
        [t,x,fs,IF_O]=signal_type_new(type);
        
        x=awgn(x,snr,'measured',lll*200+k1+3+10);
        
        
       I_max_new=FADTFD_SNR_S(x);
       %         figure;tfsapl(x,I_max_new(:,1:2:128),'grayscale','on');
        IF_COMPUTE_new_edge_link;
        %IF_COMPUTE;
        
    end
     var_snr_FADTFD(lll,:)=(mean(mse1))';
  
    for k1=1:N_S
        [t,x,fs,IF_O]=signal_type_new(type);
        
        x=awgn(x,snr,'measured',lll*200+k1+3+10);
        I_max_new=tfr_stft_high(hilbert(x'));
        % figure; imshow(1-I_max_new,[]);
        %                 figure;imshow(1-I_fixed_wind,[]);
        IF_COMPUTE_new_edge_link;
        %IF_COMPUTE;
        
    end
    var_snr_AFS(lll,:)=(mean(mse1))';
%     
%
 for k1=1:N_S
        [t,x,fs,IF_O]=signal_type_new(type);
        
        x=awgn(x,snr,'measured',lll*200+k1+3+10);
       
        %I_max_new =  S_method(x, 2,85 );
        I_max_new =  S_method(conj(x'), 15,85 );
        I_max_new(I_max_new<0)=0;
        
        
        IF_COMPUTE_new_edge_link;
      %  IF_COMPUTE;
        
    end
    var_snr_S_method(lll,:)=(mean(mse1))';
%     
   
    
    for k1=1:N_S
        [t,x,fs,IF_O]=signal_type_new(type);
        
        x=awgn(x,snr,'measured',lll*200+k1+3+10);
        wvdz=wvd(x,length(x)-1,1,2^nextpow2(length(x)));
        
             g=extnd_mbd(0.05,0.2,1,length(x));
                gsig=ifft(fft(wvdz.').');

        smg=gsig.*g;
        I_max_new=real(fft(ifft(smg.').'));
        I_max_new=imresize(I_max_new,[length(x) length(x)]);
        I_max_new(I_max_new<0)=0;
       % mesh(I_max_new)
        IF_COMPUTE_new_edge_link;
        %IF_COMPUTE;
        
        
    end
    var_snr_embd(lll,:)=(mean(mse1))';
    
     for k1=1:N_S
        [t,x,fs,IF_O]=signal_type_new(type);
        
        x=awgn(x,snr,'measured',lll*200+k1+3+10);
%         wvd=wvd(x,length(x)-1,1,2^nextpow2(length(x)));;
        [hadtfd]=HTFD_new12(x);

        hadtfd(hadtfd<0)=0;
        hadtfd=filter2(B2,hadtfd);
        I_max_new = hadtfd;

        
            I_max_new=imresize(I_max_new,[length(x) length(x)]);
        
        I_max_new(I_max_new<0)=0;
        
        IF_COMPUTE_new_edge_link;
        %IF_COMPUTE;
        
        
    end
    var_snr_HADTFD(lll,:)=(mean(mse1))';
    
   
    
    
    
%     
    for k1=1:N_S
        [t,x,fs,IF_O]=signal_type_new(type);
        
        x=awgn(x,snr,'measured',lll*200+k1+3+10);
        I_max_new=real(quadtfd(x,length(x)/2-1,1,'cw',64,128));
        I_max_new(I_max_new<0)=0;
        %I_max_new =  S_method(x, 2,85 );
%         I_max_new =  S_method(conj(x'), 15,85 );
%         I_max_new(I_max_new<0)=0;
%         
        
        IF_COMPUTE_new_edge_link;
      %  IF_COMPUTE;
        
    end
    var_snr_CW(lll,:)=(mean(mse1))';
   
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



% snr=-5:3:15;
% 
% display('frac tfd');
% log10(var_snr_AFS)
% save new_mse_AFS_type2.mat var_snr_AFS
% 
% 
% display('MBD');
% log10(var_snr_embd)
% save new_mse_mbd_type2.mat var_snr_mbd
% 
% display('HTFD');
% log10(var_snr_HTFD)
% save new_mse_mbd_type2.mat var_snr_HTFD


snr=-5:2:-1;
% for i=1:1
% %      figure;
% %     plot(snr,mean(log10(var_snr_mbd')),'-co','linewidth',3);
% % %    
%     hold on;
%     plot(snr,mean(log10(var_snr_spec')),':gs','linewidth',3);
%     
%     hold on;
%     plot(snr,mean(log10(var_snr_embd')),'-.b+','linewidth',3);
%     hold on;
%     plot(snr,mean(log10(var_snr_CSK')),'--bd','linewidth',3);
%     
%     hold on;
%     plot(snr, mean(log10(var_snr_AFS')),'-rh','linewidth',3);
% 
% %     hold on;
% %     plot(snr, mean(log10(var_snr_HTFD')),'-.b+','linewidth',3);
% 
%     xlabel('Signal to Noise Ratio');
%     ylabel('log10(Mean Square Error)');
%     axis([min(snr) max(snr)  -5  0])
% end
snr=-5:2:16;
% snr=-5:5:20;

for i=1:2
     figure;
    plot(snr,10*(log10(var_snr_AOK(:,i))),'-co','linewidth',4);
   
    hold on;
    plot(snr,10*(log10(var_snr_CW(:,i))),':gs','linewidth',4);
%     
    hold on;
    plot(snr,10*(log10(var_snr_embd(:,i))),'-.b+','linewidth',4);
   hold on;
   plot(snr,10*(log10(var_snr_HADTFD(:,i))),'--md','linewidth',4);
    
    hold on;
    plot(snr, 10*(log10(var_snr_AFS(:,i))),'-rh','linewidth',4);

    hold on;
    plot(snr, 10*(log10(var_snr_FADTFD(:,i))),'-.k+','linewidth',4);
    
    
    legend(' AOK','CW',' EMBD',' HADTFD','AFS',' FADTFD','Location','Best')
    xlabel('Signal to Noise Ratio');
    ylabel('Mean Square Error (dB)');
   % axis([min(snr) max(snr)  -50  0])
end