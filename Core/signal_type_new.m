function [ t,x,fs,IF_O ] = signal_type_new(type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


if type==1
    
    
    
    fs=256;
    t=0:1/fs:1-1/fs;
    s=1*exp(1i*(40*pi*t.^3+70*pi*t))+   0.65*exp(1i*(40*pi*t.^3+40*pi*t));
    s=real(s);
    
    IF_O(:,1)=(60*t.*t+35);
    IF_O(:,2)=(60*t.*t+20);
    IF_O=IF_O/(fs/2);
    x=s' ;
    
elseif type==2
    
    t=0:127;
    x=0.65*cos(2*pi*(0.46*t-0.000007*t.^3))+cos(2*pi*(0.4*t-0.000007*t.^3));
    IF_O(:,1)=0.46-0.000021*t.^2;
    IF_O(:,2)=0.4-0.000021*t.^2;
    IF_O=2*IF_O;
    fs=1;
    x=x';
    
elseif type==3
    fs=256;
    t=0:1/fs:1-1/fs;
    s=1*exp(1i*(40*pi*t.^3+70*pi*t))+  0.9*exp(1i*(40*pi*t.^3+50*pi*t))+0.65*exp(1i*(40*pi*t.^3+10*pi*t));
    s=real(s);
    
    IF_O(:,1)=(60*t.*t+35);
    IF_O(:,2)=(60*t.*t+25);
    IF_O(:,3)=(60*t.*t+5);
    IF_O=IF_O/(fs/2);
    x=s' ;
    
elseif type==4
    fs=1;
    t=0:1:255;
    s=1*exp(1i*(2*pi*(0.15/512)*t.^2+2*0.1*pi*t))+  1*exp(1i*(2*pi*(0.15/512)*t.^2+2*0.15*pi*t));
    s=real(s);
    
    IF_O(:,1)=(0.3*t/512+0.1);
    IF_O(:,2)=(0.3*t/512+0.15);
    IF_O=IF_O/(1/2);
    x=s' ;

  elseif type==5
    fs=1;
    t=0:1:255;
    s=1*exp(1i*2*pi(0.001*t.^2+0.1*t))+  1*exp(1i*2*pi*(0.001*t.^2+0.15*t));
    s=real(s);
    
    IF_O(:,1)=(0.002*t+0.1);
    IF_O(:,2)=(0.002*t+0.15);
    IF_O=IF_O/(1/2);
    x=s' ;

  elseif type==6
    fs=1;
    t=0:1:255;
    s=1*cos(2*pi*(0.05*t+0.000001*t.^3))+1*cos(2*pi*(0.075*t+0.000001*t.^3))+0.5*cos(2*pi*(0.45*t-0.0001*t.^2/2));
    IF_O(:,1)=(0.05+0.000003*t.^2);
    IF_O(:,2)=(0.075+0.000003*t.^2);
    IF_O(:,3)=(0.45-0.0002*t/2);
        IF_O=IF_O/(1/2);

x=s';
 elseif type==7
    fs=1;
    t=0:1:255;
    s=1*cos(2*pi*(0.025*t+0.000001*t.^3))+1*cos(2*pi*(0.05*t+0.000001*t.^3))+0.5*cos(2*pi*(0.45*t-0.00013*t.^2))+0.5*cos(2*pi*(0.475*t-0.00013*t.^2));
    IF_O(:,1)=(0.025+0.000003*t.^2);
    IF_O(:,2)=(0.05+0.000003*t.^2);
    IF_O(:,3)=0.45-0.00026*t;
        IF_O(:,4)=0.475-0.00026*t;

        IF_O=IF_O/(1/2);

x=s';
elseif type==8
    fs=1;
    t=0:1:255;
    s=1*cos(2*pi*(0.025*t+0.000001*t.^3))+1*cos(2*pi*(0.05*t+0.000001*t.^3))+0.5*cos(2*pi*(0.45*t-0.0000005*t.^3))+0.5*cos(2*pi*(0.475*t-0.0000005*t.^3));
    IF_O(:,1)=(0.025+0.000003*t.^2);
    IF_O(:,2)=(0.05+0.000003*t.^2);
    IF_O(:,3)=0.45-0.0000015*t.^2;
        IF_O(:,4)=0.475-0.0000015*t.^2;

        IF_O=IF_O/(1/2);

x=s';
elseif type==9
    fs=1;
    t=0:1:255;
    s=1*cos(2*pi*(0.025*t+0.000001*t.^3))+1*cos(2*pi*(0.05*t+0.000001*t.^3))+0.5*cos(2*pi*(0.45*t-0.0000005*t.^3))+0.5*cos(2*pi*(0.475*t-0.0000005*t.^3));%+0.5*cos(2*pi*(0.3*t));%+[zeros(1,75) cos(2*pi*(0.25*t(76:end-75))) zeros(1,75)];
    IF_O(:,1)=(0.025+0.000003*t.^2);
    IF_O(:,2)=(0.05+0.000003*t.^2);
    IF_O(:,3)=0.45-0.0000015*t.^2;
        IF_O(:,4)=0.475-0.0000015*t.^2;
        IF_O(:,5)=0.25;

        IF_O=IF_O/(1/2);

x=s';
end
