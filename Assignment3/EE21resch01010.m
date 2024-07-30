%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment : 03
% Name       : Shantanu Yadav
% Roll No    : EE20MTECH12001
% Course     : DSP Lab 2021
% 
% Details    : This file generates OFDM pulses with digital modulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%% Inputs
baseFreq = 5;             % Base Frequency
symTime = 1/baseFreq;     % Symbol Time
totalSubcarr = 6;         % Total Subcarriers (Should be less than or equal to FFT size)
fftSize = 16;             % FFT Size (Should be non zero integer)
%% 
fs=fftSize*baseFreq;
Tsubcarr=symTime/fftSize; %Time for each sub carrier
Ton=totalSubcarr*Tsubcarr;
Ts=1/fs;
t_i=0:Ts:Ts*(fftSize-1);t_i=t_i';

K=50;
T=Tsubcarr/K;
f=1/Tsubcarr;
t=0:T:symTime-T;
%% Basis function for PSK Modulation
bpsk_basis_function=sqrt(2/Ton)*cos(2*pi*f*t);
qpsk_basis_function=[sqrt(2/Ton)*cos(2*pi*f*t) ; sqrt(2/Ton)*sin(2*pi*f*t)] ;
bpsk_basis_size=size(bpsk_basis_function);
qpsk_basis_size=size(qpsk_basis_function);
%% Generate random bits with uniform distribution 
%%Map the bits to BPSK and QPSK 
% no. of the symbols equals to no. of subcarrier in assignment 2.
%generate N random numbers in the interval (a,b) with the formula r = a + (b-a).*rand(N,1). 
a=-1;
b=1;
bpsk_data=a+(b-a).*rand(totalSubcarr,1);
qpsk_data=a+(b-a).*rand(totalSubcarr,2);
x_bpsk=psk(bpsk_data,2,1)';
x_qpsk=psk(qpsk_data,4,1)';
%% Generation of OFDM pulses as in assignment 2
xk_bpsk=[zeros(1,(fftSize-totalSubcarr)/2) x_bpsk' zeros(1,(fftSize-totalSubcarr)/2)]';

xk_qpsk=[zeros(2,(fftSize-totalSubcarr)/2) x_qpsk' zeros(2,(fftSize-totalSubcarr)/2)]';

subplot(611); stem(xk_bpsk),title("BPSK Input Bits"),xlabel("n"),ylabel("xk_bpsk");
subplot(612);stem(xk_qpsk(:,1)),title("QPSK In-Phase Bits"),xlabel("n"),ylabel("xk_qpsk_inphase");
subplot(613);stem(xk_qpsk(:,2)),title("QPSK Out-Phase Bits"),xlabel("n"),ylabel("xk_qpsk_outphase");
%% Mapping PSK bits to modulated signal
%BPSK
bpsk_signal=reshape((xk_bpsk.*reshape(bpsk_basis_function,[K,fftSize])')',bpsk_basis_size);
subplot(614);plot(t,bpsk_signal);title('BPSK bits mapped to carrier'),xlabel('time(t)'),ylabel('bpsk_signal');

%QPSK
qpsk_signal_inphase=reshape((xk_qpsk(:,1).*reshape(qpsk_basis_function(1,:),[K,fftSize])')',bpsk_basis_size);
qpsk_signal_outphase=reshape((xk_qpsk(:,2).*reshape(qpsk_basis_function(2,:),[K,fftSize])')',bpsk_basis_size);
subplot(615);plot(t,qpsk_signal_inphase);title('QPSK IN-Phase bits mapped to carrier'),xlabel('time(t)'),ylabel('qpsk_signal_inphase');
subplot(616);plot(t,qpsk_signal_outphase);title('QPSK OUT-Phase bits mapped to carrier'),xlabel('time(t)'),ylabel('qpsk_signal_outphase');
%% IDFT of OFDM signal using summation
% This section of the code takes the IDFT of the OFDM signal 
si_bpsk=idft(xk_bpsk,fftSize);
si_qpsk=idft(xk_qpsk(:,1)+j*xk_qpsk(:,2),fftSize);

%Power Normalization
si_bpsk=si_bpsk/sqrt((fftSize*sum(abs(si_bpsk).^2)));
si_qpsk=si_qpsk/sqrt((fftSize*sum(abs(si_qpsk).^2)));
%% IDFT of OFDM signal using Matlab IFFT
% This section of the code takes the IFFT of the OFDM signal 

% Both the output (your IDFT and Matlab IFFT) output should exactly match
si2_bpsk=ifft(xk_bpsk,fftSize);
si2_qpsk=ifft(xk_qpsk(:,1)+j*xk_qpsk(:,2),fftSize);

%Power Normalization
si2_bpsk=si2_bpsk/sqrt((fftSize*sum(abs(si2_bpsk).^2)));
si2_qpsk=si2_qpsk/sqrt((fftSize*sum(abs(si2_qpsk).^2)));

figure(2);
subplot(411),plot(t_i,abs(si_bpsk)),title("IDFT(BPSK)"),xlabel("time(t)"),ylabel("s(t)_bpsk");
subplot(412),plot(t_i,abs(si2_bpsk)),title("IFFT(BPSK)"),xlabel("time(t)"),ylabel("s(t)_bpsk");
subplot(413),plot(t_i,abs(si_qpsk)),title("IDFT(QPSK)"),xlabel("time(t)"),ylabel("s(t)_qpsk");
subplot(414),plot(t_i,abs(si2_qpsk)),title("IFFT(QPSK)"),xlabel("time(t)"),ylabel("s(t)_qpsk");
%% Plotting the spectrum of the signal
% This section plots the spectrum of the signal

% Spectrum of the signal should exactly match with the subcarriers you have
% taken in the input. 

% If the spectrum doesn't match for any specific case then you should
% mention that case along with the reasons

% Power of the all the frequency bins should be 1. And you should mentioned
% how you have done the power adjustment (if you have done it)
Xk_bpsk=fft(si_bpsk);
Xk_qpsk=fft(si_qpsk);
Xk_qpsk=[real(Xk_qpsk) imag(Xk_qpsk)];

%Calculating power
pow_bpsk=sum(abs(Xk_bpsk).^2);
pow_qpsk=sum(abs(Xk_qpsk).^2);
disp('BPSK Power = '),disp(pow_bpsk);
disp('QPSK Power ='),disp(sum(pow_qpsk)); %in-phase power + out-phase power

%Plotting
figure(3);
subplot(311); stem(Xk_bpsk),title("Output BPSK Bits"),xlabel("n"),ylabel("Xk_bpsk");
subplot(312);stem(Xk_qpsk(:,1)),title("Output QPSK In-Phase Bits"),xlabel("n"),ylabel("Xk_qpsk_inphase");
subplot(313);stem(Xk_qpsk(:,2)),title("Output QPSK Out-Phase Bits"),xlabel("n"),ylabel("Xk_qpsk_outphase");


function X=psk(data,m,A)
%% BPSK
if (m==2)
    bpsk_const=A*[-1,1];
    X=zeros(size(data))';
    [~,I]=sort(SquareDist(bpsk_const',data));
    idx=I(1,:);
    for i=1:length(data)
        X(i)=bpsk_const(:,idx(i));
    end
end
%% QPSK
if m==4
    qpsk_const=A*[[1;1],[-1;1],[-1;-1],[1;-1]];
    X=zeros(size(data))';
    [~,I]=sort(SquareDist(qpsk_const',data));
    idx=I(1,:);
    for i=1:length(data)
        X(:,i)=qpsk_const(:,idx(i));
    end
end
end
%% Function to calculate Euclidean Distance
function D = SquareDist(X1, X2)
% SQUAREDIST - computes Euclidean SQUARED distance matrix
%
% E = distance(X1, X2)
%
%    X1 - (NxD) matrix
%    X2 - (MxD) matrix
%
% Returns:
%    E - (NxM) Euclidean SQUARED distances between vectors in X1 and X2

n = size(X1, 1);
m = size(X2, 1);

sq1 = sum(X1.*X1, 2);
sq2 = sum(X2.*X2, 2);

D = sq1*ones(1, m) + ones(n, 1)*sq2' - 2*(X1*X2');
end
%% 
function Yf=idft(yk,N)
for k=0:N-1
    sum=0;
    for j=0:N-1   
        sum=sum+yk(j+1)*exp(1i*2*pi*k*j/N);
    end
    Yf(k+1)=sum;
end
Yf=Yf.'/N;
end
