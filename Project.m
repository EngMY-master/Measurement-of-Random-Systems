

% Decimation code 

load('project.mat')


Ts = t(2); 
Fs = 1/Ts; 

%% Decimate x by a factor of 4 (BW is now 20 KHz/4 = 5 kHz

Fn_new  = Fs/4; 
% Design a low-pass Butterworth filter

[B,A] = butter(6,Fn_new/Fs); 

% compute the filter frequency response 
[H,wout] = freqz(B, A);

%freqs(B,A)                  % Plot frequency response

figure(1)
subplot(211), semilogy(wout/2/pi,abs(H))
subplot(212), plot(wout/2/pi,180/pi*angle(H))



xf = filter(B,A,x); 
xf = xf(1:4:end); 


% change the corner frequency and the skipping points 
