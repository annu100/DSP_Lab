%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment : 04
% Name       : ANNU
% Roll No    : EE21RESCH01010
% Course     : DSP Lab 2021
% 
% Details    : This file generates OFDM pulses with digital modulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%% Inputs
baseFreq = 5;             % Base Frequency
symTime = 1/baseFreq;     % Symbol Time
%totalSubcarr = 256;         % Total Subcarriers (Should be less than or equal to FFT size)
totalSubcarr=4096;
fftSize = 4096;                  % FFT Size (Should be non zero integer)
CP=128;
%cyclic prefix size
Eb=1; 
Eb_No_dB=-5:17;
SNRdB=-5:1:17;                             %Signal to Noise Ratio (in dB)
SNR=10.^(SNRdB/10);

time_offsetarray=1:400;%0 to 10 percent of fftsize
mse1=[];

estimated_timeoffset=[];
for theta=1:80:400
for count=1:length(SNR)
   N1=100; %No of OFDM symbols
   Err=0;
   for iter=1:N1
       
        %Bit generation
        databits=randi([0,1],1,2*totalSubcarr);
        %databits=randint(1,totalSubcarr);         %Generate binary data source
        %Bits mapping to bpsk 
        bits_bpsk=bpsk_map(databits);
        %Bits mapping to bpsk 
        bits_qpsk=QPSK_mod(databits);
        No=(2*Eb)/SNR(count);
        stand(count)=sqrt(No/2);   %Generate AWGN 
        
        %% QPSK
        ofdm_tx=sqrt(fftSize)*ifft(bits_qpsk,fftSize);
        %   Introduce offset
        ofdm_qpad=zeros(1,length(ofdm_tx));
        ofdm_offset_qpsk=[ofdm_qpad(:,end-theta+1:end) ofdm_tx(:,1:end-theta)];
        ofdm_qpad=ofdm_tx;
           % 5. Cyclic prefix addition qpsk
        cyclic_pre_qpsk=ofdm_offset_qpsk(:,end-CP+1:end);
        ch_in=[cyclic_pre_qpsk ofdm_offset_qpsk];
        %Noise Modelling
        noise = randn(1,length(ch_in)) + 1j*randn(1,length(ch_in));
        ch_out = ch_in + stand(count)*noise;
        %ch_out = ch_out(:,end-fftSize)
        
        
        %Remove cyclic prefix
        remov_cp=ch_out((CP+1):end);
        %Take fft
        
        
        fftout=fft(ch_out,fftSize);
        recieved=fftout;
        %for i = 1:fftSize
        %recieved(i)= exp(-1i*2*pi*baseFreq*theta*i)*recieved(i);
        %end
        x_k=fft(ofdm_tx,fftSize)
        channel=recieved./x_k;
        %phase=0;
        %for m=1:fftSize-1
        %    phase=phase+(channel(m)*conj(channel(m+1)));
        %end
        Const_array =conj(channel(:,2:end)).*channel(:,1:end-1);
        phase=sum(Const_array)/length(Const_array);
        n_est=(abs((fftSize*log(phase))/(2*pi)));
        Err = Err + (abs(n_est)-theta)^2;
        end
        %n_est(theta)=mean(n_est);
        mse1(count)=Err/N1;
end
lege = ['Offset = ',num2str(theta)];
   semilogy(SNRdB,mse1,'linewidth',1.5,'DisplayName',lege);
   hold on;
   
title("mean square error")
ylabel("SNR (in Db)")
xlabel("mean square error for different time offset")
end

legend show;
disp(n_est);
%disp(size(n_est))
%disp(size(mse1))

function data=bpsk_map(bits)
  data=2*bits-1;
end

function data=qpsk_map(bits)
  data=[];
  len=length(bits);
  for l=1:2:len
      data=[data,(1/sqrt(2))*(2*((bits(l)-0.5)+1i*(bits(l+1)-0.5)))];
  end
end
function qpsk=QPSK_mod(x)
    x_odd = x(:,1:2:end-1);
    x_even = x(:,2:2:end);
    qpsk=(1/sqrt(2))*(2*((x_odd-0.5)+1i*(x_even-0.5)));
end

