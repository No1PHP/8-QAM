clear all;
close all;
tic
M = 8;
b = log2(M);
N = 3*10^6;         
input = randi([0,1],1,N);

% reshape the input signal with 3 a group, so that a column represents
% a symbol(3 bits)
input_tri = reshape(input, [3, N/3]);

d = 1;
E_s2=(2^2-1)*d/12;
E_s4=(4^2-1)*d/12;
Es=E_s2+E_s4;

EsN0 = -2:1:30;
N0 = Es./(10.^(EsN0./10));  
EbN0_db = EsN0 - 10*log10(3);
EbN0 = 10.^(EbN0_db./10);
std_dev =reshape(sqrt(N0./2),1,1,length(EsN0)); %standard deviation

s1 = [-d/2; 3*d/2];
s2 = [d/2; 3*d/2];
s3 = [-d/2; 1*d/2];
s4 = [d/2; 1*d/2];
s5 = [-d/2; -1*d/2];
s6 = [d/2;- 1*d/2];
s7 = [-d/2; -3*d/2];
s8 = [d/2; -3*d/2];

bit1 = [0,1,0].'; 
bit2 = [1,1,0].'; 
bit3 = [0,1,1].'; 
bit4 = [1,1,1].'; 
bit5 = [0,0,1].'; 
bit6 = [1,0,1].';
bit7 = [0,0,0].'; 
bit8 = [1,0,0].'; 

s_i = [s1, s2, s3, s4, s5, s5, s6, s7, s8];  
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------ Gray mapping part
%--------------------------------------------------------------------------
% symbol mapping matrix

%--------------------------------------------------------------------------
% bit checking
b1 = repmat(bit1, 1, N/3); %map s1
b2 = repmat(bit2, 1, N/3); %map s2
b3 = repmat(bit3, 1, N/3); %map s3
b4 = repmat(bit4, 1 ,N/3); %map s4
b5 = repmat(bit5, 1, N/3); %map s5
b6 = repmat(bit6, 1, N/3); %map s6
b7 = repmat(bit7, 1, N/3); %map s7
b8 = repmat(bit8, 1, N/3); %map s8


% xor the target and the signal, return matrices 
% which have the location of the corresponding mapping signal
symbol1_mat = xor(input_tri,b1);
symbol2_mat = xor(input_tri,b2);
symbol3_mat = xor(input_tri,b3);
symbol4_mat = xor(input_tri,b4);
symbol5_mat = xor(input_tri,b5);
symbol6_mat = xor(input_tri,b6);
symbol7_mat = xor(input_tri,b7);
symbol8_mat = xor(input_tri,b8);


% logic index, to identify the mapping signals' location
% those who in the symbol_mat that equals to 0 after xor 
% is the corresponding mapping
location_s1 = ~logical(sum(symbol1_mat));
location_s2 = ~logical(sum(symbol2_mat));
location_s3 = ~logical(sum(symbol3_mat));
location_s4 = ~logical(sum(symbol4_mat));
location_s5 = ~logical(sum(symbol5_mat));
location_s6 = ~logical(sum(symbol6_mat));
location_s7 = ~logical(sum(symbol7_mat));
location_s8 = ~logical(sum(symbol8_mat));


% make it two dimensional to later point product the symbol
%i.e. generate the mapped signal

% using the location index to get the mapped signal value
signal1 = location_s1.*s1;
signal2 = location_s2.*s2;
signal3 = location_s3.*s3;
signal4 = location_s4.*s4;
signal5 = location_s5.*s5;
signal6 = location_s6.*s6;
signal7 = location_s7.*s7;
signal8 = location_s8.*s8;


signals = signal1+signal2+signal3+signal4+signal5+signal6+signal7+signal8;

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------ Add gaussian white noise part
%-------------------------------------------------------------------------
% build a 3D array, where:
% row is different index and for different noise addition, namely 
% different N0
noise =std_dev.*randn(2,length(signals),length(EsN0));

% AWGN
received = signals + noise;
% 8 signals ------ 8 decision condition
% 3 threshold 0 -d d
s1_shot = received(1,:,:)<0 & received(2,:,:)>d;
s2_shot =received(1,:,:)>0 & received(2,:,:)>d;
s3_shot = received(1,:,:)<0 & received(2,:,:)<d & received(2,:,:)>0;
s4_shot = received(1,:,:)>0 & received(2,:,:)<d & received(2,:,:)>0;
s5_shot = received(1,:,:)<0 & received(2,:,:)<0 & received(2,:,:)>-d;
s6_shot = received(1,:,:)>0 & received(2,:,:)<0 & received(2,:,:)>-d;
s7_shot = received(1,:,:)<0 & received(2,:,:)<-d;
s8_shot = received(1,:,:)>0 & received(2,:,:)<-d;

dm_s1 = s1_shot.*s1;
dm_s2 = s2_shot.*s2;
dm_s3 = s3_shot.*s3;
dm_s4 = s4_shot.*s4;
dm_s5 = s5_shot.*s5;
dm_s6 = s6_shot.*s6;
dm_s7 = s7_shot.*s7;
dm_s8 = s8_shot.*s8;

demodulated = dm_s1 + dm_s2 + dm_s3 + dm_s4 + dm_s5 + dm_s6 + dm_s7 + dm_s8;

demodulated_bits = s1_shot.*bit1 + s2_shot.*bit2 + s3_shot.*bit3 + s4_shot.*bit4 + s5_shot.*bit5 + s6_shot.*bit6 + s7_shot.*bit7 + s8_shot.*bit8;
toc
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------ Evaluation part
%-------------------------------------------------------------------------

P2 = (2*(2-1)/2).*Q(d./sqrt(2.*N0));
P4 = (2*(4-1)/4).*Q(d./sqrt(2.*N0));        
Theoretical = P2+P4-P2.*P4;

is_equal = (demodulated == signals);
symbol_error =  reshape(is_equal(1,:,:) & is_equal(2,:,:),length(received),1,length(EsN0));
symbol_error_num = squeeze(length(received) - sum(symbol_error));
figure(1);
semilogy(EsN0,symbol_error_num./(length(signals)),'mx--');

hold on 
semilogy(EsN0,Theoretical,'b.-');
xlabel('Es/N0 (dB)');
ylabel('Symbol Error Rate');
legend('SER','theoretical');
title('SER');

is_equal = reshape(demodulated_bits == input_tri,length(signals),3,length(EsN0));
bit_error_num = squeeze(N - (sum(is_equal(:,1,:))+sum(is_equal(:,2,:))+sum(is_equal(:,3,:))));
figure(2);
semilogy(EbN0_db,bit_error_num./N);
line = ones(length(EsN0),1)./10^4;
hold on
semilogy(EbN0_db,line,'m--');
xlabel('Eb/N0(dB)');
ylabel('Bit Error Rate');
title('BER');

function y=Q(x)
y=erfc(x./sqrt(2))/2;
end



