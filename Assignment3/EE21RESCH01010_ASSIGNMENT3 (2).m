%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name       : ANNU
% Roll No.   : EE21RESCH01010
% Assignment : 03
% Course     : DSP Lab 2021
% 
% Details    : This file generates OFDM pulses 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% 1.Generate random bits with uniform distribution 
clc;close all;clear all

%% Inputs
baseFreq = 5;             % Base Frequency
symTime = 1/baseFreq;     % Symbol Time
totalSubcarr = 16;         % Total Subcarriers (Should be less than or equal to FFT size)
fftSize = 16;             % FFT Size (Should be non zero integer)
%% 
fs=fftSize*baseFreq;
Ton=(totalSubcarr*symTime)/fftSize;
Ts=1/fs;

N=16;         % Number of data bits needs to be generated


% data bits generation %
data=randi([0,1],1,N);

%%%%%%% 2.Map the bits to BPSK and QPSK 
br=5; %Let us transmission bit rate is 5 hz
f=br; % minimum carrier frequency
T=1/br; % bit duration
t=0:T/100:T; % Time vector for one bit information

% X-axis for QPSK %
for i=1:N/2
    T1(:,i)=(i-1)*symTime:T/100:i*symTime;
end
% X-axis for BPSK %
for i=1:N
    T2(:,i)=(i-1)*symTime:T/100:i*symTime;
end
%%%%%%%%%% BPSK Generation %%%%%%%%%%

for i=1:N
    if data(i)==0
        psk(:,i)=-1*sin(2*pi*5*t);
        bpskData(i)=-1;
    else
        psk(:,i)=sin(2*pi*5*t);
        bpskData(i)=1;
    end
end
disp('bpsk data')
disp(psk)
%%%%%%%% Plotting BPSK %%%%%%%%
figure(1)
stem(data)
title('Data bits')
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%XXXX QPSK Modulation without consideration of noise XXXXX
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% % information
%data=randi([0,1],1,N);

data_NZR=2*data-1; % Data Represented at NZR form for QPSK modulation
s_p_data=reshape(data_NZR,2,length(data)/2);  % S/P convertion of data


br=5; %Let us transmission bit rate is 5hz
f=br; % minimum carrier frequency
T=1/br; % bit duration
t=T/99:T/99:T; % Time vector for one bit information



% XXXXXXXXXXXXXXXXXXXXXXX QPSK modulation  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
y=[];
y_in=[];
y_qd=[];
for(i=1:length(data)/2)
    y1=s_p_data(1,i)*cos(2*pi*f*t); % inphase component
    y2=s_p_data(2,i)*sin(2*pi*f*t) ;% Quadrature component
    y_in=[y_in y1]; % inphase signal vector
    y_qd=[y_qd y2]; %quadrature signal vector
    y=[y y1+y2]; % modulated signal vector
end
Tx_sig=y; % transmitting signal after modulation
tt=T/99:T/99:(T*length(data))/2;

figure(5)
%QPSK plotting 

subplot(3,1,1);
plot(tt,y_in,'linewidth',3), grid on;
title(' wave form for inphase component in QPSK modulation ');
xlabel('time(sec)');
ylabel(' amplitude');

subplot(3,1,2);
plot(tt,y_qd,'linewidth',3), grid on;
title(' wave form for Quadrature component in QPSK modulation ');
xlabel('time(sec)');
ylabel(' amplitude');


subplot(3,1,3);
plot(tt,Tx_sig,'r','linewidth',3), grid on;
title('QPSK modulated signal (sum of inphase and Quadrature phase signal)');
xlabel('time(sec)');
ylabel(' amplitude');


%%%%%%%%%% QPSK Generation %%%%%%%%%
r=1;
for i=1:2:N
     %%%%%%Time resolution
    t=0:T/100:T
    if (data(i)==0 && data(i+1)==0)
        qpsk(:,r)=sin(2*pi*5*t);
       
    elseif (data(i)==0 && data(i+1)==1)
        qpsk(:,r)=-1*sin(2*pi*5*t);
      
    elseif (data(i)==1 && data(i+1)==0)
        qpsk(:,r)=cos(2*pi*5*t);
        
    elseif (data(i)==1 && data(i+1)==1)
         qpsk(:,r)=-1*cos(2*pi*5*t);
        
    else
        disp('oops')
    end
    r=r+1;
end
disp('qpsk')
disp(qpsk)
%%%%%% BPSK and QPSK plotting%%%%%%%%
figure(2) 
subplot(2,1,1)
plot(T2,psk,'linewidth',3)
title('BPSK mapping')
subplot(2,1,2)
plot(T1,qpsk,'linewidth',3)
title('QPSK mapping(sum of inphase as well as quadrature) ')


%%%%%%%%%%SubCarriers Allocation to the data bits for BPSK%%%%%%
S2P=reshape(bpskData,N/totalSubcarr,totalSubcarr)
for i=1:size(S2P,2)
    Subcarrier(:,i)=S2P(:,i);
end

%%%%%%%%%% OFDM with BPSK %%%%%%%%%%%%%
% symbol time 'symTime'
%iZero_pad=[zeros(1,7) Subcarrier(1,2) zeros(1,7)] ;
%zero padding to generate idft effienciently
%dft=zeros(1,fftSize);
%idftdata=zeros(1,fftSize); %idft calculation
%for n=0:fftSize-1
%    for k=0:fftSize-1
%        idftdata(n+1)=idftdata(n+1)+(iZero_pad(k+1)*exp(i*2*pi*k*n/fftSize));
%    end
%end
%No zero padding is required this time as no.  of subcarriers I took is 16
idftData=idft(Subcarrier,totalSubcarr)
%%%%%%IFFT built in function%%%%%%%
ifftData=ifft(Subcarrier,N);
%%%%%Plots Comparison between Manually computed IDFT and Matlab IFFT%%%%%%%
figure(3)
subplot (2,1,1)
stem(ifftshift(abs(idftData)))
title('Manually Computed IDFT for bpsk')
subplot(2,1,2)
stem(ifftshift(abs(ifftData)))
title('Matlab IFFT for bpsk')
%%%%%%%%%% OFDM with QPSK %%%%%%%%%%%%%
%%%%%%%%%%SubCarriers Allocation to the data bits for QPSK%%%%%%
S2Pq=reshape(Tx_sig,N/totalSubcarr,[])
for i=1:size(S2Pq,2)
    Subcarrierq(:,i)=S2Pq(:,i);
end
idftDataq=idft(Subcarrierq,totalSubcarr)
%%%%%%IFFT built in function%%%%%%%
ifftDataq=ifft(Subcarrierq,N);
%%%%%Plots Comparison between Manually computed IDFT and Matlab IFFT%%%%%%%
figure(6)
subplot (2,1,1)
stem(ifftshift(abs(idftDataq)))
title('Manually Computed IDFT for QPSK')
subplot(2,1,2)
stem(ifftshift(abs(ifftDataq)))
title('Matlab IFFT for qpsk')
%%%%%%%Spectrum of the ofdm Pulses%%%%%
dftData=fft(ifftData,totalSubcarr);
figure(4)
subplot(2,1,1)
stem(bpskData)
title('actual data bits in Subcarriers i.e BPSK')
subplot(2,1,2)
stem(dftData)
title('Spectrum of OFDM Data ')
dftDataq=fft(ifftDataq,totalSubcarr);
figure(7)
subplot(2,1,1)
stem(Tx_sig)
title('actual data bits in Subcarriers i.e QPSK')
subplot(2,1,2)
stem(dftDataq)
title('Spectrum of OFDM Data using QPSK(both inphase and quadrature)')
%Power normalisation
pow_bpsk=sum((abs(bpskData).^2)/16);
pow_qpsk=sum((abs(Tx_sig).^2)/ 792);
disp('BPSK Power = '),disp(pow_bpsk);
disp('QPSK Power ='),disp(sum(pow_qpsk)); %in-phase power + out-phase power

%%%%% bit Energies for all bits are Normalized to 1 %%%%%%
%%%%% all bits carriers unit Energy %%%%%
%%%%% Power of all Freq bins is 1 as amplitude of bits is 1 or 0 %%%%% 
%%%%% No power adjustment is done since amplitudes of bits are either 1 or zero  %%%%%
%%%%% FFT shift is needed to properly view the data points and it is giving
%%%%% appropriate results also.

function Yf=idft(yk,N) %function for idft calculation
for k=0:N-1
    sum=0;
   for i=0:N-1   
      sum=sum+yk(i+1)*exp(j*2*pi*k*i/N);
   end
    Yf(k+1)=sum;
end
Yf=Yf'/N;
end

% XXXXXXXXXXXXXXXXXXXXXXXXX    end of program    XXXXXXXXXXXXXXXXXXXXXXXXXX
    
