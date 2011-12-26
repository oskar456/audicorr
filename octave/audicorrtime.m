clear all;
close all;

samplerate = 8000;
needle = wavread('needle.wav');
haystack = wavread('haystack.wav');

n1 = needle(:,1);
h1 = haystack(:,1);

nenergy = sum(n1.^2); % Normalize to needle's energy

len = length(h1)-length(n1)+1;

c = zeros(len, 1);
for i=1:len
    for j=1:length(needle)
        c(i) =  c(i) + n1(j) * h1(j+i-1);
    end
end
        
c = c ./ nenergy;
timescale = (0:(len-1))./samplerate;
plot(timescale,c);

[m, i] = max(c);

fprintf(1, 'Maximalni korelace v hodnote %d je v case %d sekund\n', m, timescale(i));
