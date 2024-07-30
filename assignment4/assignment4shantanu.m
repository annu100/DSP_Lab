%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment : 04
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
%totalSubcarr = 256;  % Total Subcarriers (Should be less than or equal to FFT size)
totalSubcarr=2^10;
%totalSubcarr=16;
fftSize = 256;     % FFT Size (Should be non zero integer)
%fftSize = 4;                  % FFT Size (Should be non zero integer)
%L=32                       %cyclic prefix size
L=32;
%L=2;
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

bpsk_bits=randi([0,1],totalSubcarr,1);
qpsk_bits=randi([0,1],totalSubcarr,2);

bpsk_data=2*bpsk_bits-1;
qpsk_data=(2*qpsk_bits(:,1)-1)+1i*(2*qpsk_bits(:,2)-1);

%% OFDM Transmitter
parallel_bpsk=reshape(bpsk_data,[fftSize,totalSubcarr/fftSize ]);
parallel_si_bpsk=idft_mat(parallel_bpsk,fftSize);
si_bpsk=reshape(parallel_si_bpsk,[totalSubcarr,1]);
si_bpsk_cp=[si_bpsk(totalSubcarr-L+1:totalSubcarr)' si_bpsk']';

parallel_qpsk=reshape(qpsk_data,[fftSize,totalSubcarr/fftSize ]);
parallel_si_qpsk=idft_mat(parallel_qpsk,fftSize);
si_qpsk=reshape(parallel_si_qpsk,[totalSubcarr,1]);
si_qpsk_cp=[si_qpsk(totalSubcarr-L+1:totalSubcarr)' si_qpsk']';
%% AWGN noise
%channel model: y=x+n where n is complex noise
M_bpsk=length(si_bpsk_cp);
M_qpsk=length(si_qpsk_cp);

Es_bpsk=sum(abs(si_bpsk_cp).^2)/M_bpsk;
%Es_bpsk=sum(abs(bpsk_data).^2)/M_bpsk;
Eb_bpsk=Es_bpsk/1;
Es_qpsk=sum(abs(si_qpsk_cp).^2)/M_qpsk;
Eb_qpsk=Es_qpsk/2;

SNR=1
for i=0:10
    SNRb(i+1)=i*SNR;
    SNRs_bpsk(i+1)=SNRb(i+1)+10*log(1);
    SNRs_qpsk(i+1)=SNRb(i+1)+10*log(2);
    
    SNR_lin_bpsk(i+1)=10^(SNRs_bpsk(i+1)/10);
    SNR_lin_qpsk(i+1)=10^(SNRs_qpsk(i+1)/10);
    
    sigma_bpsk(i+1)=sqrt(Es_bpsk/(2*SNR_lin_bpsk(i+1)));
    sigma_qpsk(i+1)=sqrt(Es_qpsk/(2*SNR_lin_qpsk(i+1)));
    
    n_bpsk=sigma_bpsk(i+1)*(randn(M_bpsk,1)+1i*randn(M_bpsk,1));
    n_qpsk=sigma_qpsk(i+1)*(randn(M_qpsk,1)+1i*randn(M_qpsk,1));
    
    r_bpsk_cp=si_bpsk_cp+n_bpsk;
    r_qpsk_cp=si_qpsk_cp+n_qpsk;
    
%% OFDM Receiver
    r_bpsk=r_bpsk_cp(L+1:M_bpsk,:);
    r_qpsk=r_qpsk_cp(L+1:M_qpsk,:);
    
    parallel_r_bpsk=reshape(r_bpsk,[fftSize,totalSubcarr/fftSize ]);
    parallel_r_qpsk=reshape(r_qpsk,[fftSize,totalSubcarr/fftSize ]);
    
    parallel_X_bpsk=fft(parallel_r_bpsk,fftSize);
    parallel_X_qpsk=fft(parallel_r_qpsk,fftSize);
    
    X_bpsk=reshape(parallel_X_bpsk,[totalSubcarr,1]);
    X_qpsk=reshape(parallel_X_qpsk,[totalSubcarr,1]);
%% ML Detector
    estimated_bpsk=psk(X_bpsk,2,1)';
    estimated_bpsk_bits=(estimated_bpsk+1)/2;
    
    estimated_qpsk=psk([real(X_qpsk) imag(X_qpsk)],4,1)';
    estimated_qpsk_bits=(estimated_qpsk+[1 1])/2;
    
%% BER
    BER_bpsk(i+1)=sum(bpsk_bits~=estimated_bpsk_bits)/totalSubcarr;
    theoretical_BER_bpsk(i+1)=qfunc(sqrt(2*SNR_lin_bpsk(i+1)));
    
    err_qpsk=([real(qpsk_data) imag(qpsk_data)]~=estimated_qpsk);
    BER_qpsk(i+1)=(sum(err_qpsk(:,1))+sum(err_qpsk(:,2)))/(2*totalSubcarr);
    theoretical_BER_qpsk(i+1)=1*qfunc(sqrt(SNR_lin_qpsk(i+1)));
end
figure(1);
hold on
semilogy(SNRb,BER_bpsk,'b-*');
semilogy(SNRb,theoretical_BER_bpsk,'r-');
%axis([-1 9 -0.5 1])
%axis([0 9 .99e-4 1]);
xlabel('SNR (dB)');
ylabel('BER');
legend('BER bpsk','Theoretical bpsk');

figure(2);
hold on
semilogy(SNRb,BER_qpsk,'b-*');
semilogy(SNRb,theoretical_BER_qpsk,'r-');
legend('BER qpsk','Theoretical qpsk');
%axis([-1 9 .99e-4 1]);





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
%% IDFT using Matrix
function xn=idft_mat(Xk,N)
%Generating NxN idft mat
IDFT_mat=zeros(N,N);
W=exp(1i*2*pi/N);
for n=0:N-1
    for k=0:N-1
        IDFT_mat(n+1,k+1)=(W^(n*k))/N;
    end
end
%display(IDFT_mat);
xn=IDFT_mat*Xk;
end
