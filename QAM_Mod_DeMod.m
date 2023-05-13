clc;
clear;
close all;
% Generate binary data
N = 10000; % Length of data
data = randi([0 1], N, 1); % Random binary data stream
% Modulate using 16-QAM
M = 16; % Number of symbols in 16-QAM
qam16 = qammod(data, M, 'InputType', 'bit', 'UnitAveragePower', true); % Modulated signal
% Plot constellation diagram
figure;
plot(real(qam16), imag(qam16), 'o');
xlabel('In-phase'); ylabel('Quadrature');
title('16-QAM constellation diagram');
% Demodulate QAM signal back to binary data
demodulated_data = qamdemod(qam16, M, 'OutputType', 'bit', 'UnitAveragePower', true);
% Compute bit error rate (BER)
num_errors = sum(xor(data, demodulated_data));
ber = num_errors/N;
disp(['Bit error rate (BER): ', num2str(ber)]);
