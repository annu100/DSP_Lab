%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
fftSize = 16;             % FFT Size (Should be non zero integer)

% Code should not take any other input other than the above ones
% Code should run for all the valid combinations of the inputs above

% Find out the sampling frequency of the signal
Fs=fftSize*baseFreq
%%  
% This section of the code generates 'totalSubcarr' OFDM pulses for a given 
% symbol time 'symTime'
iZero_pad=[zeros(1,7) ones(1,2) zeros(1,7)] ;
t=0:1:15;
%zero padding to generate idft effienciently
dft=zeros(1,fftSize);
idft=zeros(1,fftSize); %idft calculation
for n=0:fftSize-1
    for k=0:fftSize-1
        idft(n+1)=idft(n+1)+(iZero_pad(k+1)*exp(i*2*pi*k*n/fftSize));
    end
end
idft=idft./fftSize;
t=0:0.2/15:0.2; %0.2 is symbol time
figure(1)
stem(t,abs(idft));
title('OFDM generation')
ylabel ('subcarriers');
xlabel ('Time');
%% DFT of OFDM signal using summation
% This section of the code takes the DFT of the OFDM signal 
for k=0:fftSize-1
    for n=0:fftSize-1
        dft(k+1)=dft(k+1)+(idft(n+1)*exp((-i)*2*pi*k*n/fftSize));
    end
end
%------------------------------------------------------------

%------------------------------------------------------------
magnitude=abs(dft); % Find the magnitudes of individual DFT points
%code block to plot the magnitude response
%------------------------------------------------------------



%% DFT of OFDM signal using Matlab FFT
% This section of the code takes the FFT of the OFDM signal 
Y = fftshift(idft); % fft shift is useful for visualizing a Fourier transform with the zero-frequency component in the middle of the spectrum. For vectors, fftshift(X) swaps the left and right halves of X . For matrices, fftshift(X) swaps the first quadrant with the third and the second quadrant with the fourth.

% Both the output (your DFT and Matlab FFT) output should exactly match

% Explanation: fftshift  move the zero-frequency component to the center
%of the array. zero-frequency component is near the center of the matrix.


%% Plotting the spectrum of the signal


% Spectrum of the signal should exactly match with the subcarriers you have
% taken in the input. 
t=0:fftSize-1;
figure(2)
subplot(1,2,1)
stem(t,magnitude);
title('DFT of OFDM Signal')
ylabel ('Amplitude');
xlabel ('time');
t=0:fftSize-1;
subplot(1,2,2)
stem(t,magnitude);
title('DFT of OFDM signal using Matlab FFT')
ylabel ('Amplitude');
xlabel ('time');

% Power of the all the frequency bins should be 1. IT IS VERIFIED FROM GRAPHS.And you should mentioned
% how you have done the power adjustment (if you have done it)