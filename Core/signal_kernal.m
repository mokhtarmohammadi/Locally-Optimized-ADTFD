function analytic_sig_ker = signal_kernal(x);
% Matlab programs to develop points learned in the understanding of
% time-frequency distributions
%
% This old girl is a signal kernel function.
%    Both real and analytic versions.
%
%   K(n,m) = z(n+m)z*(n-m)
%
% where, z() is the analytic associate of real input signal s()
%    n is the sampled or discrete time
%    m is the sample or disrete lag
%
%  Example functions and their outputs
%
%  REAL
%
%
%
%
%
%  ANALYTIC
%
%
% Nathan Stevenson
% April 2004 
%

N = length(x);
%real_x = x;

if mod(length(x),2) == 0
    true_X = fft(x);
    analytic_X = [true_X(1) 2.*true_X(2:N/2) true_X(N/2+1) zeros(1,N/2-1)];
    analytic_x = ifft(analytic_X);    
else
    true_X = fft(x);
    analytic_X = [true_X(1) 2.*true_X(2:ceil(N/2)) zeros(1,floor(N/2))];    
    analytic_x = ifft(analytic_X);    
end

analytic_sig_ker = zeros(N,N);
for m = -round(N/2-1):1:round(N/2-1);
    analytic_sig_ker(m+round(N/2)+1,:) = sig_ker_corr(analytic_x,m); 
%    real_sig_ker(m+round(N/2)+1,:) = sig_ker_corr(real_x,m);
end

