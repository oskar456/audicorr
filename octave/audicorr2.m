clear all;
close all;

samplerate = 8000;
needle = wavread('needle.wav');
haystack = wavread('haystack.wav');

n1 = needle(:,1);
h1 = haystack(:,1);

nenergy = sum(n1.^2); % Normalize to needle's energy

len = max(length(n1), length(h1));
nfft = 2^nextpow2(len);
spech = fft(h1,nfft);
spech = spech(1:nfft/2+1);
specn = fft(n1,nfft);
specn = specn(1:nfft/2+1);
convspec = spech .* conj(specn);
convspec = [convspec;; conj(convspec(end-1:-1:2)) ]; %Restore hermitian symmetry
c = ifft( convspec );
% We only search for positive lag
c = c(1:len)./nenergy;

timescale = (0:(len-1))./samplerate;
plot(timescale,c);

[m, i] = max(c);

fprintf(1, 'Maximalni korelace v hodnote %d je v case %d sekund\n', m, timescale(i));