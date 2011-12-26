clear all;
close all;

samplerate = 8000;
needle = wavread('needle.wav');
haystack = wavread('haystack.wav');

n1 = needle(:,1);
h1 = haystack(:,1);

nenergy = sum(n1.^2); % Normalize to needle's energy

len = max(length(n1), length(h1));

c = xcorr(h1, n1);
% We only search for positive lag
c = c(len:len*2-1)./nenergy;

timescale = (0:(len-1))./samplerate;
plot(timescale,c);

[m, i] = max(c);

fprintf(1, 'Maximalni korelace v hodnote %d je v case %d sekund\n', m, timescale(i));