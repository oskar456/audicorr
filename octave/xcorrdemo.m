
x = rand(100,1);
len = length(x);

%# autocorrelation
nfft = 2^nextpow2(2*len-1);
r = ifft( fft(x,nfft) .* conj(fft(x,nfft)) );

%# rearrange and keep values corresponding to lags: -(len-1):+(len-1)
r = [r(end-len+2:end) ; r(1:len)];

%# compare with MATLAB's XCORR output
all( (xcorr(x)-r) < 1e-10 )
