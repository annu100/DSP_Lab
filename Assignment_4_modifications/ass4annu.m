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

total_com=2^20; %I request sir to keep 2^27 which which implies 524288 iterations.
%iterations = 2^20/256=4096;  % Total Subcarriers (Should be less than or equal to FFT size)
fftSize = 256;     % FFT Size (Should be non zero integer)
                 % FFT Size (Should be non zero integer)
%L=32                       %cyclic prefix size
Cyclic_len=32;

%% 
fs=fftSize*baseFreq;
Tsubcarr=symTime/fftSize; %Time for each sub carrier

Ts=1/fs;
Eb_No_dB=0:10;
for Eb_No_dB=0:10
    n=1;    %Es=nEb where n=1 for BPSK
    Es_No_dB=Eb_No_dB+10*log10(n);  
    Es_No=10^( Es_No_dB/10);
    
    BPSK_symbols=[1,-1];    %BPSK constellation set (Already Normalised)
    E_sb=1/2*sum(BPSK_symbols.^2);
    sigmab = sqrt(E_sb/(2*Es_No));
    
end
fprintf('Normalised bpsk energy')
disp(E_sb)
    %QPSK
for Eb_No_dB=0:10
    n=2;    %Es=nEb where n=2 for QPSK
    Es_No_dB=Eb_No_dB+10*log10(n);  
    Es_No=10^( Es_No_dB/10);
    
    QPSK_symbols=(1/sqrt(2))*[1+1i,1-1i,-1+1i,-1-1i];    
    E_sq=1/4*sum(abs(QPSK_symbols).^2); %Normalising QPSK constellation set for unit energy
    sigmaq = sqrt(E_sq/(2*Es_No));
end
fprintf('Normalised qpsk energy')
disp(E_sq)
%%No iterations for ofdm symbol required,as all ofdm symbols,computations
%%are being done together in order to minimise time and computation 
%%This is done using matrix form using matlab ifft and fft 
bdata_bits=randi([0,1],total_com,1);
qdata_bits=randi([0,1],total_com,2);
%BPSK MAPPING
bpsk_data=bpsk_map(bdata_bits);
qpsk_data=(2*qdata_bits(:,1)-1)+1i*(2*qdata_bits(:,2)-1);
%qpsk_data=qpsk_map(bpsk_bits)
%% OFDM Transmitter
%Serial to parallel conversion for bpsk
s2p_tx_bpsk=reshape(bpsk_data,[fftSize,total_com/fftSize ]);
%IDFT USING MATRIX
ofdm_si_bpsk=ifft(s2p_tx_bpsk,fftSize);
%PARALLEL TO SERIAL CONVERSION for bpsk
p2s_tx_bpsk=reshape(ofdm_si_bpsk,[total_com,1]);
%CYCLIC PREFIX ADDITION
cp_bpsk=[p2s_tx_bpsk(total_com-Cyclic_len+1:total_com)' p2s_tx_bpsk']';
%SERIAL TO PARALLEL CONVERSION for qpsk
s2p_tx_qpsk=reshape(qpsk_data,[fftSize,total_com/fftSize ]);
%IDFT USING MATRIX
ofdm_si_qpsk=ifft(s2p_tx_qpsk,fftSize);
%Parallel to serial
p2s_tx_qpsk=reshape(ofdm_si_qpsk,[total_com,1]);
%Cyclic prefix addition
cp_qpsk=[p2s_tx_qpsk(total_com-Cyclic_len+1:total_com)' p2s_tx_qpsk']';
%% AWGN noise
%channel model: y=x+n where n is complex noise
lenbpsk=(length(cp_bpsk));
lenqpsk=(length(cp_qpsk));

Es_bpsk=sum(abs(cp_bpsk).^2)/(length(cp_bpsk));
%Es_bpsk=sum(abs(bpsk_data).^2)/M_bpsk;
Eb_bpsk=Es_bpsk/1; %Since n=1 for bpsk
Es_qpsk=sum(abs(cp_qpsk).^2)/(length(cp_qpsk));
Eb_qpsk=Es_qpsk/2;

SNR=1
for i=0:10
    %standard deviation for qpsk and bpsk is findout separately calculated
    SNRb(i+1)=i*SNR;
    SNRs_bpsk(i+1)=SNRb(i+1)+10*log(1);
    SNRs_qpsk(i+1)=SNRb(i+1)+10*log(2);
    
    SNR_lin_bpsk(i+1)=10^(SNRs_bpsk(i+1)/10);
    SNR_lin_qpsk(i+1)=10^(SNRs_qpsk(i+1)/10);
    
    sigma_bpsk(i+1)=sqrt(Es_bpsk/(2*SNR_lin_bpsk(i+1)));
    sigma_qpsk(i+1)=sqrt(Es_qpsk/(2*SNR_lin_qpsk(i+1)));
    %NOISE IS GAUSSIAN ADDITIVE NOISE
    n_bpsk=sigma_bpsk(i+1)*(randn((length(cp_bpsk)),1)+1i*randn((length(cp_bpsk)),1));
    n_qpsk=sigma_qpsk(i+1)*(randn((length(cp_qpsk)),1)+1i*randn((length(cp_qpsk)),1));
    
    r_bpsk_cp=cp_bpsk+n_bpsk;
    r_qpsk_cp=cp_qpsk+n_qpsk;
    
%% OFDM Receiver
%CYCLIC PREFIX REMOVAL BEFORE IFFT
    r_bpsk=r_bpsk_cp(Cyclic_len+1:(length(cp_bpsk)),:);
    r_qpsk=r_qpsk_cp(Cyclic_len+1:(length(cp_qpsk)),:);
    %Serial to parallel
    s2p_r_bpsk=reshape(r_bpsk,[fftSize,total_com/fftSize ]);
    %Parallel to serial
    p2s_r_qpsk=reshape(r_qpsk,[fftSize,total_com/fftSize ]);
    %fft
    s2p_X_bpsk=fft(s2p_r_bpsk,fftSize);
    %fft
    p2s_X_qpsk=fft(p2s_r_qpsk,fftSize);
    %parallel to serial
    databits_bpsk=reshape(s2p_X_bpsk,[total_com,1]);
    %parallel to serial
    databits_qpsk=reshape(p2s_X_qpsk,[total_com,1]);
%% detection using ML criteria
%ESTIMATATION AND THEN BER CALCULATION
    est_bpsk=psk(databits_bpsk,2,1)';
    est_bpsk_bits=(est_bpsk+1)/2;
    
    est_qpsk=psk([real(databits_qpsk) imag(databits_qpsk)],4,1)';
    %% 
    estimated_qpsk_bits=(est_qpsk+[1 1])/2;
    
%% BER
    BER_bpsk(i+1)=sum(bdata_bits~=est_bpsk_bits)/total_com;
    th_BER_bpsk(i+1)=qfunc(sqrt(2*SNR_lin_bpsk(i+1)));
    
    err_qpsk=([real(qpsk_data) imag(qpsk_data)]~=est_qpsk);
    BER_qpsk(i+1)=(sum(err_qpsk(:,1))+sum(err_qpsk(:,2)))/(2*total_com);
    th_BER_qpsk(i+1)=1*qfunc(sqrt(SNR_lin_qpsk(i+1)));
end
subplot(1,2,1);
hold on
semilogy(SNRb,BER_bpsk,'k-*');
semilogy(SNRb,th_BER_bpsk,'r-*');
xlabel('SNR (dB)');
ylabel('BER');
legend('BER bpsk','Theoretical bpsk');

subplot(1,2,2);
hold on
semilogy(SNRb,BER_qpsk,'k-*');
semilogy(SNRb,th_BER_qpsk,'r-*');
legend('BER qpsk','Theoretical qpsk');

%BPSK MAPPING
function data=bpsk_map(bits)
  data=2*bits-1;
end
%QPSK MAPPING
function data=qpsk_map(bits)
  data=[];
  len=length(bits);
  for l=1:2:len
      data=[data,(1/sqrt(2))*(2*((bits(l)-0.5)+1i*(bits(l+1)-0.5)))];
  end
end




%FUNCTION FOR MODULATION OF PSK AND QPSK
function mod_data=psk(databits,mod,Amp)
%% BPSK
if (mod==2)
    b_const=Amp*[-1,1];
    mod_data=zeros(size(databits))';
    [~,I]=sort(eu_dist(b_const',databits));
    idx=I(1,:);
    for i=1:length(databits)
        mod_data(i)=b_const(:,idx(i));
    end
end
%% QPSK
if (mod==4)
    q_const=Amp*[[1;1],[-1;1],[-1;-1],[1;-1]];
    mod_data=zeros(size(databits))';
    [~,I]=sort(eu_dist(q_const',databits));
    idx=I(1,:);
    for i=1:length(databits)
        mod_data(:,i)=q_const(:,idx(i));
    end
end
end
%% Function to calculate Euclidean Distance
function Dist = eu_dist(X1, X2)
p = size(X1, 1);
q = size(X2, 1);

d1 = sum(X1.*X1, 2);
d2 = sum(X2.*X2, 2);

Dist = d1*ones(1, q) + ones(p, 1)*d2' - 2*(X1*X2');
end