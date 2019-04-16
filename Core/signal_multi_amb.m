function [x,diff_x1,diff_x2]=signal_multi_amb(t,popup_amb_val)
%   SIGNAL_MULTI - MULICOMPONENT TEST SIGNALS FOR AMBIGUITY DOMAIN
%   REALIZATION  
%   [x, diff_x1, diff_x2 ]= SIGNAL_MULTI  (t,varargin)
%
%   Input parameters:  t - time domain
%                      varargin - coefficients of the signal in the case of
%                      creating new signal
%   Output parameters: x -  signal
%                      diff_x1 - the first derivative of the first part of signal x
%                      diff_x2 - the first derivative of the second part of signal x
%   Paper:"A Virtual Instrument for Time-Frequency Analysis of Signals with
%   Highly Non-Stationary Instantaneous Frequency"
%   Authors: Irena Orovic, Milica Orlandic, Srdjan Stankovic, Zdravko
%   Uskokovic
%   Date: March 2010.
% Copyrigth (c) by authors


if (popup_amb_val==1)
    x=1*exp(j*1*(4*cos(pi*t)+2/3*cos(6*pi*t)+1/2*cos(1*pi*t)+8.5*pi*t))+1*exp(j*2*(1*cos(pi*t)+1/2*cos(4*pi*t)+1/4*cos(2*pi*t))-9.5*j*pi*t);
    diff_x1=-4*pi*sin(pi*t)-2/3*6*pi*sin(6*pi*t)-1/2*pi*sin(pi*t)+7.5*pi;
    diff_x2=-2*pi*sin(4*pi*t)-1/2*4*pi*sin(4*pi*t)-1/4*2*pi*sin(2*pi*t)-9.5*pi;
elseif (popup_amb_val==2)
    x=exp(j*1.2*(3*cos(pi*t)+2/3*cos(7*pi*t)+8.5*pi*t))+1*exp(j*2*(1*cos(pi*t)+1/4*cos(2*pi*t)+1/2*cos(6*pi*t)-6.5*pi*t));
    diff_x1=1.2*(-3*pi*sin(pi*t)-2/3*7*pi*sin(7*pi*t)+8.5*pi);
    diff_x2=2*(-pi*sin(pi*t)-1/4*2*sin(2*pi*t)-1/2*6*sin(6*pi*t)-6.5*pi);
elseif (popup_amb_val==3)
fs = 16;
n = 0:1/fs:8-1/fs; N = length(n);
sig1 = 2*cos(2*pi*n*0.05*fs);
sig2 = 2*cos(2*pi*n*0.1*fs);
sigm = 0.0020;
B = fir1(16*2,0.08,'high');
s1 = 10*exp(-(fs^2*(n-0.9375).^2)/sigm);
s1 = filter(B,1,s1);
s2 = 10*exp(-(fs^2*(n-2.1875).^2)/sigm);%.*cos(0.1*pi*n);
s2 = filter(B,1,s2);
s3 = 10*exp(-(fs^2*(n- 3.4375).^2)/sigm);%.*cos(0.1*pi*n);
s3 = filter(B,1,s3);
s4 = 10*exp(-(fs^2*(n-4.6875).^2)/sigm);%.*cos(0.1*pi*n);
s4 = filter(B,1,s4);
s5 = 10*exp(-(fs^2*(n- 5.9375).^2)/sigm);%.*cos(0.1*pi*n);
s5 = filter(B,1,s5);
s = s1 + s2 + s3 + s4 + s5;
Signal = sig1 + sig2 + s;
x=Signal;
elseif (popup_amb_val==4)
    s1=2*cos(0.1*pi*t+0.0008*pi*t.^2);
s2=2*cos(0.18*pi*t+0.0008*pi*t.^2);
s3=2*cos(0.15*pi*t+pi*t.^2);
s4=2*cos(0.2*pi*t+pi*t.^2);
x=s1+s2+s3+s4;
end