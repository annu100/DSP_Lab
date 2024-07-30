 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name       : Annu
% Roll No.   : EE21RESCH01010
% Assignment : 02
% Course     : DSP Lab 2021
% 
% Details    : This file generates OFDM pulses 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%% Inputs:
baseFreq = 5;             % Base Frequency
symTime = 1/baseFreq;     % Symbol Time
totalSubcarr = 2;         % Total Subcarriers (Should be less than or equal to FFT size)
fftSize= 16;             % FFT Size (Should be non zero integer)

% Code should not take any other input other than the above ones
% Code should run for all the valid combinations of the inputs above
% Find out the sampling frequency of the signal

sampling_freq=fftSize*baseFreq

%% Generation of OFDM pulses
% This section of the code generates 'totalSubcarr' OFDM pulses for a given 
% symbol time 'symTime'
X=ones(totalSubcarr,1) % Input assumption
%sub_carriers=fftSize.*ifft(X,fftSize)
x=ones(16,1)
sub_carriers=compute_ifft(x)
% Since ofdm pulses are(IFFT *fftSize) of modulated symbols
%IFFT of y can be calculated as 
%y(n+1) = x(n+1).*exp((1j*2*pi*k*n)/N);
% where 1/N is the frequency and using k*(1/N) ,we are getting different
% frequency
figure(1)
stem(sub_carriers)
ylabel("Subcarriers of OFDM")
title("OFDM GENERATION")

%% DFT of OFDM signal using summation
% This section of the code takes the DFT of the OFDM signal 

dft_ofdm=compute_dft(sub_carriers)
fprintf("dft of ofdm pulses is %c",dft_ofdm)

%% DFT of OFDM signal using Matlab FFT
% This section of the code takes the FFT of the OFDM signal 
fft_ofdm=fft(sub_carriers)
fft_ofdm=fftshift(fft_ofdm)
%Y = fftshift( X ) rearranges a Fourier transform X by shifting the zero-frequency component to the center of the array.
%So,we can easilty compare by frequency by frequency as already shown in
%the same grapgh
fprintf("dft of ofdm pulses is %c",fft_ofdm)

% Both the output (your DFT and Matlab FFT) output should exactly match
% You should explain why you have taken fftshift and what happens if you
% don't take it
fprintf("from the plot as well as values,we can easily see that both the fft and dft of ofdm subcarriers matches and giving the same result")

%% Plotting the spectrum of the signal
% This section plots the spectrum of the signal

% Spectrum of the signal should exactly match with the subcarriers you have
% taken in the input. 
figure(2)
subplot(2,1,1)
stem(abs(dft_ofdm))
title("dft of ofdm ")
ylabel("dft of ofdm using ofdm formulla")
subplot(2,1,2)
stem(abs(dft_ofdm))
title("fft of ofdm using fft formula")
% If the spectrum doesn't match for any specific case then you should
% mention that case along with the reasons
%ANS-spectrum matches 
% Power of the all the frequency bins should be 1. And you should mentioned
% how you have done the power adjustment (if you have done it)
N=16;
Power_dft=(abs(dft_ofdm).^2)/(N.^2)	%Power of the DFT
Power_fft=(abs(fft_ofdm).^2)/(N.^2)%Power of the FFT

function [xdft]=compute_dft(x) %dft computation function 
N=length(x);
y=zeros(1,N);
xdft=zeros(1,N);
for k=0:N-1
   for n = 0:N-1
     y(n+1) = x(n+1).*exp(-(1j*2*pi*k*n)/N);
   end
   xdft(k+1)= sum(y);
end
end
function [xdft]=compute_ifft(x) %ifft computation function %OFDM PULSE GENERATION
symTime=1/5;
N=length(x);
y=zeros(1,N);
xdft=zeros(1,N);
for k=0:N-1
   for n = 0:N-1
     y(n+1) = x(n+1).*exp((1j*2*pi*k*n)*(1/symTime));
   end
   xdft(k+1)= sum(y);
end
end



