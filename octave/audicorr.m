clear all;
close all;

samplerate = 8000;
needle = wavread('needle.wav');
haystack = wavread('haystack.wav');

n1 = needle(:,1);
h1 = haystack(:,1);

nenergy = sum(n1.^2); % Normalize to needle's energy

len = max(length(n1), length(h1));
nfft = 2^nextpow2(2*len-1);
c = ifft( fft(h1,nfft) .* conj(fft(n1,nfft)) );
% We only search for positive lag
c = c(1:len)./nenergy;

timescale = (0:(len-1))./samplerate;
plot(timescale,c);

[m, i] = max(c);

fprintf(1, 'Maximalni korelace v hodnote %d je v case %d sekund\n', m, timescale(i));