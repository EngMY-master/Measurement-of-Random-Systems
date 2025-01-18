y = xx ; 
nfft = length(y); 
fsamp = 1000; 
wndw = 0; 
novlap = 0; 
% Make sure that input is a column vector
argc = size(y) ;
if (argc(1)==1); y = y' ; end
% Calculate number of available ensembles
npts = length(y) ;
ensembles = floor((npts-nfft)/(nfft-novlap))+1 ;
% Initialize ensemble indexing variables
n1 = 1 ;
n2 = nfft ;
dn = nfft-novlap ;
% Initialize psd summation storage variable
arg_sum = zeros([nfft 1]) ;
% Main program loop
for k=1:ensembles
    arg_y = y(n1:n2) ; % Extract current ensemble points
    arg_y = arg_y - mean(arg_y) ; % Remove mean
    if (wndw ~= 0)
        arg_y = arg_y.*hann(nfft) ; % Apply window if required
    end
    arg_fft = fft(arg_y,nfft) ; % FFT of ensemble
    arg_abs = abs(arg_fft).^2 ; % Modulus squared of FFT
    %arg_abs = abs(conj(arg_fft).*arg_fft) ;
    arg_sum = arg_sum + arg_abs ; % Accumulate in summation variable
    n1 = n1+dn ; % Increment ensemble index variables
    n2 = n2+dn ;
end % End of main loop
arg_sum = arg_sum/ensembles ; % Average value of summed spectra
% Compute window function normalization factor
if (wndw ~= 0)
    wndw_sc = sum(hann(nfft).^2) ;
else
    wndw_sc = nfft ;
end

% Power spectrum is symmetric about Nyquist frequency, use lower half and
% multiply by 4. Then divide by 2 to convert peak^2 to rms^2. Thus scale
% factor is 2.
p = 2*arg_sum(1:nfft/2+1) ;
% Normalize to correct for window
p = p/wndw_sc ;
% Normalize spectral density for sampling units
p = p/fsamp ;
% Create frequency vector
df = fsamp/nfft ;
f = [0:df:fsamp/2]';
% Calculate overall rms level from area under PSD curve
oarms = sqrt(sum(p.*df)) ;