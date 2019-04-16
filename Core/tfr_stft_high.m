function [tfd] = tfr_stft_high(s)

L=15;
L=10;
M=20;
time=2;
II_g=zeros([time*length(s) length(s) L]);
II_hn=II_g;
II_hm=II_g;
II_rect=II_g;
II_hn1=II_g;
II_hm1=II_g;
II_rect1=II_g;

i=0;
%w1=gausswin(41,2);
for k=0:1:L-1
    i=i+1;
    w=gausswin(i*12+1+M,2);
    %w=gausswin(i*58+1+M,2);
    %w=fracft(w1,1*k/L);
    %w=ones(1,i*14+1+30-30)';
    w_1=fracft(w, 0.05);%0.05
    w_2=fracft(w, 0.1);%0.1
    w_3=fracft(w, 0.15);%0.15

    
   [ II_g(:,:,i)]=tfrsp(conj(s)',1:length(s),time*length(s), w);%/sc;
   [ II_hn(:,:,i)]=tfrsp(conj(s)',1:length(s),time*length(s),w_1);%/sc;
   [ II_hm(:,:,i)]=tfrsp(conj(s)',1:length(s),time*length(s),w_2);%/sc;
  [II_rect(:,:,i)]=tfrsp(conj(s)',1:length(s),time*length(s),w_3);%/sc;

   
%    II_g(:,:,i)=stft1(conj(s(1:600))',200,600, w);%/sc;
%    II_hn(:,:,i)=stft1(conj(s(1:600))',200,600,w_1);%/sc;
%    II_hm(:,:,i)=stft1(conj(s(1:600))',200,600,w_2);%/sc;
%    II_rect(:,:,i)=stft1(conj(s(1:600))',200,600,w_3);%/sc;
    w_1=fracft(w, -0.05);%0.05
    w_2=fracft(w, -0.1);%0.1
    w_3=fracft(w, -0.15);%0.15

   [ II_hn1(:,:,i)]=tfrsp(conj(s)',1:length(s),time*length(s),w_1);%/sc;
   [ II_hm1(:,:,i)]=tfrsp(conj(s)',1:length(s),time*length(s),w_2);%/sc;
   [ II_rect1(:,:,i)]=tfrsp(conj(s)',1:length(s),time*length(s),w_3);%/sc;

%   II_g(:,:,i)=stft1(conj(s(1:600))',600, w);%/sc;
%    II_hn1(:,:,i)=stft1(conj(s(1:600))',200,600,w_1);%/sc;
%    II_hm1(:,:,i)=stft1(conj(s(1:600))',200,600,w_2);%/sc;
%    II_rect1(:,:,i)=stft1(conj(s(1:600))',200,600,w_3);%/sc;
%    
end  
k=0;
for j=1:1:L
    k=k+1;
    for i=1:length(s)
%     II_g(:,i,k)=II_g(:,i,k)/sum(II_g(:,i,k));
%     II_hn(:,i,k)=II_hn(:,i,k)/sum(II_hn(:,i,k));
%     II_hm(:,i,k)=II_hm(:,i,k)/sum(II_hm(:,i,k));
%     II_rect(:,i,k)=II_g(:,i,k)/sum(II_rect(:,i,k));
%     II_hn1(:,i,k)=II_hn1(:,i,k)/sum(II_hn1(:,i,k));
%     II_hm1(:,i,k)=II_hm1(:,i,k)/sum(II_hm1(:,i,k));
%     II_rect1(:,i,k)=II_rect1(:,i,k)/sum(II_rect1(:,i,k));
% % % %     
    
    
    end
end

  I_max1=max(max(max(max(II_g,[],3),max(II_hn,[],3)),max(II_hm,[],3)),max(II_rect,[],3));
I_max2=max(max(max(max(II_g,[],3),max(II_hn1,[],3)),max(II_hm1,[],3)),max(II_rect1,[],3));
tfd1=max(I_max1,I_max2);
%tfd1=II_g(:,:,5);
tfd=tfd1(1:end/2,:);
%tfd=tfd1;
% k=0;
% for j=1:1:L
%     k=k+1;
%     for i=1:length(s)
%     II_g(:,i,k)=II_g(:,i,k)/sum(II_g(:,i,k));
%     II_hn(:,i,k)=II_hn(:,i,k)/sum(II_hn(:,i,k));
%     II_hm(:,i,k)=II_hm(:,i,k)/sum(II_hm(:,i,k));
%     II_rect(:,i,k)=II_g(:,i,k)/sum(II_rect(:,i,k));
%     II_hn1(:,i,k)=II_hn1(:,i,k)/sum(II_hn1(:,i,k));
%     II_hm1(:,i,k)=II_hm1(:,i,k)/sum(II_hm1(:,i,k));
%     II_rect1(:,i,k)=II_rect1(:,i,k)/sum(II_rect1(:,i,k));
%     
%     
%     
%     end
% end
% % % 
%   I_max1=min(min(min(min(II_g,[],3),min(II_hn,[],3)),min(II_hm,[],3)),min(II_rect,[],3));
% I_max2=min(min(min(min(II_g,[],3),min(II_hn1,[],3)),min(II_hm1,[],3)),min(II_rect1,[],3));
% tfd=min(I_max1,I_max2);%.*tfd;

% tfd=max(II_g,[],3);
% tfd1=zeros(size(tfd));
% tfd1(20:end-20,20:end-20)=tfd(20:end-20,20:end-20);
% tfd=tfd1;