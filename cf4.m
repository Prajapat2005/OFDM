clear all

nFFT = 64; % fft size
nDSC = 52; % number of data subcarriers
nbps = 52; % number of bits per OFDM symbol (same as the number of subcarriers for BPSK)
nSym = 10^4; % number of symbols
snrdB = [0:35]; % bit to noise ratio
SNR = snrdB + 10*log10(nDSC/nFFT) + 10*log10(64/80); % converting to symbol to noise ratio

for ii = 1:length(snrdB)

    % Transmitter
    ipBit = rand(1, nbps * nSym) > 0.5; % random 1's and 0's
    ipMod = 2 * ipBit - 1; % BPSK modulation 0 --> -1, 1 --> +1
    ipMod = reshape(ipMod, nbps, nSym).'; % grouping into multiple symbols

    % Assigning modulated symbols to subcarriers from [-26 to -1, +1 to +26]
    xF = [zeros(nSym, 6) ipMod(:, [1:nbps/2]) zeros(nSym, 1) ipMod(:, [nbps/2+1:nbps]) zeros(nSym, 5)];

    % Taking FFT, the term (nFFT/sqrt(nDSC)) is for normalizing the power of transmit symbol to 1
    xt = (nFFT/sqrt(nDSC)) * ifft(fftshift(xF.')).';

    % Appending cyclic prefix
    xt = [xt(:, [49:64]) xt];

    % Multipath channel
    nTap = 10;
    ht = 1/sqrt(2) * 1/sqrt(nTap) * (randn(nSym, nTap) + j * randn(nSym, nTap));

    % Computing and storing the frequency response of the channel, for use at receiver
    hF = fftshift(fft(ht, 64, 2));

    % Convolution of each symbol with the random channel
    for jj = 1:nSym
        xht(jj, :) = conv(ht(jj, :), xt(jj, :));
    end
    xt = xht;

    % Concatenating multiple symbols to form a long vector
    xt = reshape(xt.', 1, nSym * (80 + nTap - 1));

    % Gaussian noise of unit variance, 0 mean
    nt = 1/sqrt(2) * [randn(1, nSym * (80 + nTap - 1)) + j * randn(1, nSym * (80 + nTap - 1))];

    % Adding noise, the term sqrt(80/64) is to account for the wasted
    yt = sqrt(80/64)*xt + 10^(-SNR(ii)/20)*nt;

% Receiver
yt = reshape(yt.', 80+nTap-1, nSym).'; % formatting the received vector into symbols
yt = yt(:, [17:80]); % removing cyclic prefix

% Converting to frequency domain
yF = (sqrt(nDSC)/nFFT)*fftshift(fft(yt.')).';

% Equalization by the known channel frequency response
yF = yF./hF;

% Extracting the required data subcarriers
yMod = yF(:, [6+[1:nbps/2] 7+[nbps/2+1:nbps] ]);

% BPSK demodulation
% +ve value --> 1, -ve value --> -1
ipModHat = 2*floor(real(yMod/2)) + 1;
ipModHat(find(ipModHat > 1)) = +1;
ipModHat(find(ipModHat < -1)) = -1;

% Converting modulated values into bits
ipBitHat = (ipModHat + 1) / 2;
ipBitHat = reshape(ipBitHat.', nbps * nSym, 1).';

% Counting the errors
nErr(ii) = size(find(ipBitHat - ipBit), 2);

end

simBer = nErr/(nSym*nbps);
EbNOLin = 10.^(snrdB/10);
theoryBer = 0.5 * (1 - sqrt(EbNOLin./(EbNOLin + 1)));

close all; figure
semilogy(snrdB, theoryBer, 'bs-', 'LineWidth', 2);
hold on
semilogy(snrdB, simBer, 'mx-', 'LineWidth', 2);
axis([0 35 10^-5 1])
grid on

legend('Rayleigh-Theory', 'Rayleigh-Simulation');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('BER for BPSK using OFDM in a 10-tap Rayleigh channel')