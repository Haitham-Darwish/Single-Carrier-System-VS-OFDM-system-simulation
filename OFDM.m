clc;clear;%close all;
NumberofBits = 64000;
Data=randi([0,1],[NumberofBits, 1]);
No = 10.^[-10:0.1:3];
%======================================================================%
% For Quadri-Phase Shift Keying (QPSK)
num_rep=3;
Cp=16;% Cyclic extension
FFT_points=64;
Data_rep=repelem(Data,num_rep);
Data_rep=[Data_rep; zeros(ceil(length(Data_rep)/(FFT_points*4))*FFT_points*4-length(Data_rep),1)];

% Interleaver
r_qpsk=8;
c_qpsk=16;
for i=1:r_qpsk*c_qpsk:length(Data_rep)
    temp=reshape(Data_rep(i:i+r_qpsk*c_qpsk-1), r_qpsk, c_qpsk)';
    Data_rep(i:i+r_qpsk*c_qpsk-1)=temp(:);
end


i=1;
SI=zeros(length(Data_rep)/2, 1);
SQ=zeros(length(Data_rep)/2, 1);
da=zeros(length(Data_rep)/2, 2);
% Mapper
for d=1:2:length(Data_rep),
    da(i, :)=Data_rep(d:d+1);

    if da(i, 1)==0,
        SI(i)=-1;
    else
        SI(i)=1;
    end

    if da(i, 2)==0,
        SQ(i)=-1;
    else
        SQ(i)=1;
    end

    i=i+1;
end
S_QPSK=SI+SQ*j;

% 64-IFFT block
inv_fft_qpsk=ifft(reshape(S_QPSK, [FFT_points, size(S_QPSK,1)/FFT_points]),FFT_points);
% Cyclic extensions
for i=1:size(inv_fft_qpsk, 2)
    cyclic_ex_qpsk(:,i)=[inv_fft_qpsk(end-Cp+1:end,i); inv_fft_qpsk(:,i)];
end
% channel
Noise=sqrt(permute(repmat(No/2,1,1,size(cyclic_ex_qpsk,2)),[1 3 2])).* ...
    (randn([size(cyclic_ex_qpsk), length(No)])+randn([size(cyclic_ex_qpsk), length(No)])*j);
v1=randn(size(cyclic_ex_qpsk));%+randn(size(cyclic_ex_qpsk))*j;
v2=randn(size(cyclic_ex_qpsk));%+randn(size(cyclic_ex_qpsk))*j;
R=sqrt(v1.^2+v2.^2)/sqrt(2);
%R=normalize(sqrt(v1.^2+v2.^2)/sqrt(2));
X_QPSK=R.*(cyclic_ex_qpsk)+Noise;
%% QPSK RX
% Equalizer 
QPSK_RX=(X_QPSK./R);
% Removing cyclic extension
for ii=1:size(QPSK_RX,3)
    for i=1:size(QPSK_RX, 2)
        QPSK_RX1(:,i,ii)=QPSK_RX(Cp+1:end, i,ii);
        %[inv_fft_qpsk(end-Cp+1:end,i); inv_fft_qpsk(:,i)];
    end
end
%QPSK_RX1=QPSK_RX(Cp+1:end, :);
% FFT
QPSK_RX2=fft(QPSK_RX1, FFT_points);
QPSK_RX2=reshape(QPSK_RX2, [size(S_QPSK,1),size(QPSK_RX2,3)]);
QPSK_RX=zeros(length(Data_rep), length(No));
QPSK_RX(1:2:length(Data_rep), :)=real(QPSK_RX2)>0 ;
QPSK_RX(2:2:length(Data_rep), :)=imag(QPSK_RX2)>0;

% de-inverlever
for no=1:length(No)
    for i=1:r_qpsk*c_qpsk:size(QPSK_RX,1)
        temp=reshape(QPSK_RX(i:i+r_qpsk*c_qpsk-1, no), c_qpsk, r_qpsk)';
        QPSK_RX(i:i+r_qpsk*c_qpsk-1, no)=temp(:);
    end
end

% remove repeations if any
if num_rep>1
    X=zeros(size(QPSK_RX, 1)/num_rep,size(QPSK_RX, 2), num_rep);
    s=size(QPSK_RX);
    for k=1:s(2)
      jj=1;
      for i=1:num_rep:s(1)
        X(jj, k, :)=QPSK_RX(i:i+num_rep-1, k);
        jj=jj+1;
      end
    end
    QPSK_RX=mode(X,3);
end

QPSK_BER=sum(QPSK_RX~=Data)/NumberofBits;

% Get average symbola dand bit energy
avg_symbol_energy=4*(1+1)/4;
Eb_QPSK=avg_symbol_energy/2;
theoritical_error_QPSK=1/2*erfc(sqrt(Eb_QPSK./No));
hold on
semilogy(10*log10(Eb_QPSK./No),QPSK_BER, 'r');
%semilogy(10*log10(Eb_QPSK./No),theoritical_error_QPSK, 'r');
hold off
title('BER vs EB/No');
xlabel('Eb/No in dB'); ylabel('BER');


% For Quadrature Amplitude Modulation (QAM)

Data_rep=repelem(Data,num_rep);
Data_rep=[Data_rep; zeros(ceil(length(Data_rep)/FFT_points)*FFT_points-length(Data_rep),1)];
% Interleaver
r_qam=16;
c_qam=16;
for i=1:r_qam*c_qam:length(Data_rep)
    temp=reshape(Data_rep(i:i+r_qam*c_qam-1), r_qam, c_qam)';
    Data_rep(i:i+r_qam*c_qam-1)=temp(:);
end


i=1;
SI_QAM=zeros(length(Data_rep)/4, 1);
SQ_QAM=zeros(length(Data_rep)/4, 1);
da=zeros(length(Data_rep)/4, 4);

for d=1:4:length(Data_rep)
    da(i, :)=Data_rep(d:d+3);

    if da(i, 1)==0 && da(i, 2)==0,
        SI_QAM(i)=-3;
    elseif da(i, 1)==0 && da(i, 2)==1,
        SI_QAM(i)=-1;
    elseif da(i, 1)==1 && da(i, 2)==1,
        SI_QAM(i)=1;
    elseif da(i, 1)==1 && da(i, 2)==0,
        SI_QAM(i)=3;
    end

    if da(i, 3)==0 && da(i, 4)==0,
        SQ_QAM(i)=-3;
    elseif da(i, 3)==0 && da(i, 4)==1,
        SQ_QAM(i)=-1;
    elseif da(i, 3)==1 && da(i, 4)==1,
        SQ_QAM(i)=1;
    elseif da(i, 3)==1 && da(i, 4)==0,
        SQ_QAM(i)=3;
    end



    i=i+1;
end
S_QAM=SI_QAM+SQ_QAM*j;

% 64-IFFT block
inv_fft_qam=ifft(reshape(S_QAM, [FFT_points, size(S_QAM,1)/FFT_points]),FFT_points);
% Cyclic extensions
for i=1:size(inv_fft_qam, 2)
    cyclic_ex_qam(:,i)=[inv_fft_qam(end-Cp+1:end,i); inv_fft_qam(:,i)];
end
% channel
Noise=sqrt(permute(repmat(No/2,1,1,size(cyclic_ex_qam,2)),[1 3 2])).* ...
    (randn([size(cyclic_ex_qam), length(No)])+randn([size(cyclic_ex_qam), length(No)])*j);
v1=randn(size(cyclic_ex_qam))+randn(size(cyclic_ex_qam))*j;
v2=randn(size(cyclic_ex_qam))+randn(size(cyclic_ex_qam))*j;
%R=sqrt(v1.^2+v2.^2)/sqrt(2);
R=normalize(sqrt(v1.^2+v2.^2)/sqrt(2));
X_QAM=R.*(cyclic_ex_qam)+Noise;
%% QPSK RX
% Equalizer 
QAM_RX=(X_QAM./R);
% Removing cyclic extension
for ii=1:size(QAM_RX,3)
    for i=1:size(QAM_RX, 2)
        QAM_RX1(:,i,ii)=QAM_RX(Cp+1:end, i,ii);
        %[inv_fft_qpsk(end-Cp+1:end,i); inv_fft_qpsk(:,i)];
    end
end
%QPSK_RX1=QPSK_RX(Cp+1:end, :);
% FFT
QAM_RX2=fft(QAM_RX1, FFT_points);
QAM_RX2=reshape(QAM_RX2, [size(S_QAM,1),size(QAM_RX2,3)]);

avg_symbol_energy=(4*(3^2+3^2)+4*(1+1)+8*(3^2+1))/16;
Eb_QAM=avg_symbol_energy/4;
theoritical_error_QAM=3/8*erfc(sqrt(Eb_QAM./(2.5*No)));

QAM_RX=zeros(length(Data_rep), length(No));
QAM_RX(1:4:length(Data_rep), :)=real(QAM_RX2)>0;
QAM_RX(2:4:length(Data_rep), :)=abs(real(QAM_RX2))<2;

QAM_RX(3:4:length(Data_rep), :)=imag(QAM_RX2)>0;
QAM_RX(4:4:length(Data_rep), :)=abs(imag(QAM_RX2))<2;

for no=1:length(No)
    for i=1:r_qam*c_qam:size(QAM_RX,1)
        temp2=reshape(QAM_RX(i:i+r_qam*c_qam-1, no), c_qam, r_qam)';
        QAM_RX(i:i+r_qam*c_qam-1, no)=temp2(:);
    end
end

if num_rep>1
    Y=zeros(size(QAM_RX, 1)/num_rep,size(QAM_RX, 2), num_rep);
    q=size(QAM_RX);
    for k2=1:q(2)
      jj=1;
      for ii=1:num_rep:q(1)
        Y(jj, k2, :)=QAM_RX(ii:ii+num_rep-1, k2);
        jj=jj+1;
      end
    end
    QAM_RX=mode(Y,3);
end




QAM_BER=sum(QAM_RX~=Data)/NumberofBits;

hold on
semilogy(10*log10(Eb_QAM./No),QAM_BER, 'k');
%semilogy(10*log10(Eb_QAM./No),theoritical_error_QAM, 'k');
hold off
title('BER vs EB/No');
xlabel('Eb/No in dB'); ylabel('BER');



legend('BER of QPSK fading no rep', 'BER of QAM fading no rep', 'BER QPSK with rep', 'QAM BER with rep');
%legend('BER of QAM fading', 'BER QPSK AWGN', 'BER of QAM fading', 'QAM BER AWGN');

title('BER vs EB/No');
xlabel('Eb/No in dB'); ylabel('BER');
%yticks(10.^[-6:1:0])
xticks([0:2:24])
%ylim([10^-1 10^0])
xlim([0 18])
%set(gca, 'YScale', 'log')