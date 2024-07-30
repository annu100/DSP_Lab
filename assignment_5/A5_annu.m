%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment : 05
% Name       : ANNU
% Roll No    : EE21RESCH01010
% Course     : DSP Lab 2021
% 
% Details    : TIME OFFSET ESTIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;
%%This program is for checking using different approach i.e ML estimation
%%of offset using cyclic prefix instead of pilots
%% Parameter declaration according to requirements
fftSize = 512;
CP = 32;
Nsym = 100;
baseFreq =5;
%As mentioned by sir in class theta=cp/2,perfect estimation is obtained and
%this program is for verfication only.
theta = CP/2; %for checking only 
Eb=1; 

estimated_timeoffset=[];

SNRdB=-5:20 ;                           %Signal to Noise Ratio (in dB)
SNR=10.^(SNRdB/10); %In linear



for count=1:length(SNR)

%%OFDM SYMBOL GENERATION-QPSK
nBitPerSym=2;
carriersQ=1:nBitPerSym;
ipBit=rand(1,Nsym*fftSize)>0.5; %random 1 and 0's.
oddData=ipBit(1:2:end);
evenData=ipBit(2:2:end);
qpskModulated=sqrt(1/2)*(1i*(2*oddData-1)+(2*evenData-1)); %QPSK MAPPING
Tx_q = zeros(1,Nsym*(fftSize+CP));
OFDMsym_q = zeros(1,fftSize);  
for sym = 1:Nsym/2
    OFDMsym_q = ifft(qpskModulated(fftSize*(sym-1)+1:(fftSize*sym)),fftSize)*sqrt(fftSize);
    Tx_q((fftSize+CP)*(sym-1)+1:(fftSize+CP)*sym) = [OFDMsym_q(fftSize-CP+1:fftSize) OFDMsym_q];
end

%% AWGN channel

snrq = (2*Eb)/SNR(count);
standq(count)=sqrt(snrq/2);

%%Recieved signal after adding noise.-QPSK
noiseq = standq(count)*(randn(1,Nsym*(fftSize+CP))+1i*randn(1,Nsym*(fftSize+CP)));
Rx_q = exp(1i*2*pi*baseFreq*(0:length(Tx_q)-1)./fftSize).*Tx_q + noiseq; 
%%I have estimated time offset using by exploting the use of cyclic prefix.



%%  estimation of timing offset-QPSK
theta=CP/2;
PHI_sumq = zeros(1,Nsym*(fftSize+CP)-fftSize);
GM_sumq = zeros(1,Nsym*(fftSize+CP)-fftSize);
for n = theta:Nsym*(fftSize+CP)-(fftSize+CP)
    PHIq=0;GMq=0;
    for m = n:n+CP-1    
        PHIq = PHIq+ (Rx_q(m)*conj(Rx_q(m)) + Rx_q(m+fftSize)*conj(Rx_q(m+fftSize)));
        GMq = GMq+ Rx_q(m)*conj(Rx_q(m+fftSize));    
    end
    PHI_sumq(n) = abs(GMq)- (snrq/(snrq+1))*PHIq;
    GM_sumq(n) = -angle(GMq)/(2*pi);


estimated_timeoffsetq(count)= mean(PHI_sumq);
end
end

for count=1:length(SNR)
mse1(count)=(mean((estimated_timeoffsetq(count)-theta).^2))/fftSize;


end


%% Estimation results display
figure(1);

plot(PHI_sumq);title('Estimation of timing offset-QPSK');grid on;
%%Mean square error grapgh-BPSK
figure(2);
semilogy(SNRdB,mse1,'b')

ylabel('mean square error for time offset-using BPSK')
title("rms error on time offset ")
xlabel('SNR')
legend('actual offset=CP/2')
lgd.FontSize=0.5;
