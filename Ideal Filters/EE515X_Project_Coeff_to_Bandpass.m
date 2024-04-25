clear, clc, close all
%% Butterworth
% g=[1 2 1 1];
% g=[.7654 1.8478 1.8478 .7654 1];
% g=[.61809 1.6180 2 1.6180 .6180 1]
%% Chebychev
%g=[1.5963 1.0967 1.5963 1];
%g=[1.6703 1.1926 2.3661 .8419 1.9841];
%g=[1.7058 1.2296 2.5408 1.2296 1.7058 1];
%% Maxamaly Flat
g=[1.255 .5528 .1922 1];
%g=[1.0598 .5116 .3181 .1104 1];
%g=[.9303 .4577 .3312 .209 .0718 1];
%% User Variabler
f=10e9;
bw=500e6;
R0=50;%ohms
%% Quiuck Calc
bww=2*pi*bw;
w0=2*pi*f;
%% Impedence & Frequency Scaling with a Conversion to Bandpass
%pg 414 3ed Pozar
delta=((w0+bw)-(w0-bw))/w0;

for i=1:length(g)-1
    if rem(i,2)==1
        L(i)=g(i)*R0/(w0*delta);
        C(i)=delta/(w0*g(i)*R0);
    else
        L(i)=delta*R0/(w0*g(i));
        C(i)=g(i)/(delta*w0*R0);
    end
end
disp(L/1e-9)
disp(C/1e-12)










