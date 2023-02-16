clc;clear;%close all;
NumberofBits = 36000;
Data=randi([0,1],[NumberofBits, 1]);
No = 10.^[-10:0.1:3];
%======================================================================%
% For Quadri-Phase Shift Keying (QPSK)
num_rep=1;

Data_rep=repelem(Data,num_rep);

% Interleaver
for i=1:16:length(Data_rep)
    temp=reshape(Data_rep(i:i+15), 4, 4)';
    Data_rep(i:i+15)=temp(:);
end


i=1;
SI=zeros(length(Data_rep)/2, 1);
SQ=zeros(length(Data_rep)/2, 1);
da=zeros(length(Data_rep)/2, 2);
Noise=sqrt(No/2).*randn(length(Data_rep)/2, length(No))+sqrt(No/2).*randn(length(Data_rep)/2, length(No))*j;
v1=randn(length(Data_rep)/2, 1);%+randn(length(Data_rep)/2, 1)*j;
v2=randn(length(Data_rep)/2, 1);%+randn(length(Data_rep)/2, 1)*j;
R=sqrt(v1.^2+v2.^2)/sqrt(2);
%R=normalize(sqrt(v1.^2+v2.^2)/sqrt(2));
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
X_QPSK=R.*(S_QPSK)+Noise;




QPSK_RX=zeros(length(Data_rep), length(No));
QPSK_RX(1:2:length(Data_rep), :)=real(X_QPSK.*(conj(R)./norm(R)))>0 ;
QPSK_RX(2:2:length(Data_rep), :)=imag(X_QPSK.*(conj(R)./norm(R)))>0;


for no=1:length(No)
    for i=1:16:size(QPSK_RX,1)
        temp=reshape(QPSK_RX(i:i+15, no), 4, 4)';
        QPSK_RX(i:i+15, no)=temp(:);
    end
end

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
semilogy(10*log10(Eb_QPSK./No),QPSK_BER, 'r*');
semilogy(10*log10(Eb_QPSK./No),theoritical_error_QPSK, 'r');
hold off
title('BER vs EB/No');
xlabel('Eb/No in dB'); ylabel('BER');

% For Quadrature Amplitude Modulation (QAM)
i=1;
SI_QAM=zeros(length(Data_rep)/4, 1);
SQ_QAM=zeros(length(Data_rep)/4, 1);
da=zeros(length(Data_rep)/4, 4);
Noise=sqrt(No/2).*randn(length(Data_rep)/4, length(No))+sqrt(No/2).*randn(length(Data_rep)/4, length(No))*j;
v1=randn(length(Data_rep)/4, 1);%+randn(length(Data_rep)/4, 1)*j;
v2=randn(length(Data_rep)/4, 1);%+randn(length(Data_rep)/4, 1)*j;
R=sqrt(v1.^2+v2.^2)/sqrt(2);
%R=normalize(sqrt(v1.^2+v2.^2)/sqrt(2));

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
X_QAM=R.*S_QAM+Noise;

avg_symbol_energy=(4*(3^2+3^2)+4*(1+1)+8*(3^2+1))/16;
Eb_QAM=avg_symbol_energy/4;
theoritical_error_QAM=3/8*erfc(sqrt(Eb_QAM./(2.5*No)));

QAM_RX=zeros(length(Data_rep), length(No));
QAM_RX(1:4:length(Data_rep), :)=real(X_QAM./R)>0;
QAM_RX(2:4:length(Data_rep), :)=abs(real(X_QAM./R))<2;

QAM_RX(3:4:length(Data_rep), :)=imag(X_QAM./R)>0;
QAM_RX(4:4:length(Data_rep), :)=abs(imag(X_QAM./R))<2;

for no=1:length(No)
    for i=1:16:size(QAM_RX,1)
        temp2=reshape(QAM_RX(i:i+15, no), 4, 4)';
        QAM_RX(i:i+15, no)=temp2(:);
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
semilogy(10*log10(Eb_QAM./No),QAM_BER, 'ko');
semilogy(10*log10(Eb_QAM./No),theoritical_error_QAM, 'k');
hold off
title('BER vs EB/No');
xlabel('Eb/No in dB'); ylabel('BER');




legend('BER of QPSK fading', 'BER QPSK AWGN', 'BER of QAM fading', 'QAM BER AWGN');

title('BER vs EB/No');
xlabel('Eb/No in dB'); ylabel('BER');
yticks(10.^[-6:1:0])
xticks([0:2:24])
ylim([10^-6 10^0])
xlim([0 18])
set(gca, 'YScale', 'log')
