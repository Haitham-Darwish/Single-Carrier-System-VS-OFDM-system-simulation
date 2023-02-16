# Single-Carrier-System-VS-OFDM-system-simulation
The purpose of the project is to get introduced to the simulation of the single and  multicarrier communication systems
The requirements of the project are described in the following sections. 

# 1. Single Carrier System
![image](https://user-images.githubusercontent.com/68920161/219426943-28617e53-8a69-44b8-bc8a-b62a6018d34a.png)


## 1.1 Coding
Two cases are considered, no coding or rate 1/3 repetition code.
## 1.2 Interleaver
The Interleaver size is 4 by 4.
## 1.3 The Mapper 
The mapper takes the I/P data bits and produces the symbols to be transmitted on the 
channel. The modulation schemes under consideration are the QPSK and the 16QAM 
systems. Figure 2 shows the constellations.
## 1.4 The channel
The channel that will be simulated is the flat Rayleigh fading channel. For this channel 
model, the received signals y(n) is given by
y(n)=R(n)x(n)+v(n)
where x(n) is the transmitted signal, v(n) is the AWGN, and R(n) is the Rayleigh fading 
envelope. R(n) can be generated using the equation

![image](https://user-images.githubusercontent.com/68920161/219428669-50b29091-ea5d-4fd1-9139-1e7120b5a37f.png)

Where v(n) is AWGN. Note that R(n) in this case has power=1.

![image](https://user-images.githubusercontent.com/68920161/219428029-f1d66773-f2f0-4633-8606-7c05b67380bc.png)

## 1.5 The receiver
The simple receiver in the model under consideration will take the output of the channel 
and decide on the symbol transmitted. The output bit stream of the receiver is compared to 
the input bit stream and the BER is calculated. 
## 1.6 Mandatory Tasks
It is required to plot curves for the BER Vs Eb/No. Note that for the fading channel, Eb/No 
will be average Eb/No, as the fading magnitude will vary from one sample to the other. 
The other requirement is to evaluate the performance using repetition code. This is done 
by transmitting every “1” as three “1’s” and every “0” as three “0’s”. Draw BER curves in 
case of Rayleigh fading. Two figures are needed, one for QPSK (coding and no coding) 
and one for 16 QAM (coding and no coding).

# 2 OFDM system simulation
![image](https://user-images.githubusercontent.com/68920161/219428816-e516421d-1199-4555-9565-d5ac5dabe3c2.png)

## 2.1 Coding
Two cases are considered, no coding or rate 1/3 repetition code. Note that you have to 
adjust the number of input bits per OFDM symbol when using repetition code. For example 
if you use QPSK, only 21 data bits will be used per OFDM symbol, a zero is added to the 
encoded data to have 64 bits at the input of the mapper before the IFFT block.
## 2.2 Interleaver
For QPSK, the size of the interleaver is 8 by 16.
For 16QAM, the interleaver size is 16 by 16. 
## 2.3 Mapper
The mappers used are the same as those used in the single carrier system in section 1.1
## 2.4 IFFT
Use a size 64 IFFT block. In Matlab use the command “ifft”
## 2.5 Cyclic Extension 
16 samples cyclic prefix is to be added. 
## 2.6 Channel
Two channel models should be considered
a- AWGN channel: Same as single carrier system 
b- Frequency selective Fading channel: assume a 2-path fading channel h=[0.4 0 0.26 
0 0 0.4 0 0.6 0 0.5];
Receiver
Design a receiver to receive the signal described above in the two cases of AWGN and 
fading channels. Assume perfect channel knowledge at receiver.
## 2.7 Requirements
Same as in the single carrier system for the two channel models and for the coding/no-coding scenarios. Four figures will be needed. 
## 2.8 Water-filling 
Consider an OFDM system with 16 subcarriers. The signal is to be transmitted through the 
channel h=[0.4 0 0.26 0 0 0.4 0 0.6 0 0.5];
Consider a system with SNR gap =2 and the noise per subcarrier=1.
The total power available is 200, use Matlab to calculate the power allocation per 
subcarrier to maximize the transmission rate. 
The deliverable should be the full Matlab code and the intermediate and final results of 
running the code.
