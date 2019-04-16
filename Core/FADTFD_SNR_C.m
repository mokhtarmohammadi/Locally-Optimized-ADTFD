function I_max_new = FADTFD_SNR_C(x)
I_max_new = HTFD_new1(x, 2,30,82);
I_max_new(I_max_new<0)=0;
I_max_new=imresize(I_max_new,[length(x) length(x)]);
