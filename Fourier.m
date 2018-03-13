function [f,A] = Fourier(F,y)
Fs = F;                    % Sampling frequency
                  % y is signal to be analysed. 
L=max(size(y));                   % Length of signal
           % Time vector
% Sum of a 50 Hz sinusoid and a 120 Hz sinusoid 
 
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
plot(f,2*abs(Y(1:floor(NFFT/2+1)))) 

xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
A=2*abs(Y(1:floor(NFFT/2+1)));
