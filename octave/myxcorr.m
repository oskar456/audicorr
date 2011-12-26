function r = myxcorr(x,y)
len = max(length(x), length(y));
nfft = 2^nextpow2(2*len-1);
r = ifft( fft(x,nfft) .* conj(fft(y,nfft)) );

%# rearrange and keep values corresponding to lags: -(len-1):+(len-1)
r = [r(end-len+2:end) ; r(1:len)];
