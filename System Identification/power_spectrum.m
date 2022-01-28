%% Power Spectrum from FFT and Welch's Method
%% Create a Signal
% Signal parameters
srate = 1000; % sampling rate in Hz
time  = -2:1/srate:2;
npnts = length(time);
frex  = [ 5 8 16 20 27 ]; 
ampl  = [ 7 6 2  4  1 ]; 

% Initialize signal
signal = zeros(1,length(time));

for fi=1:length(frex)
    AM = ampl(fi)*interp1(randn(20,1),linspace(1,20,length(time)),'pchip');
    signal = signal + AM.*sin(2*pi*frex(fi)*time);
end

% Let's look at how the signal looks like
figure(1)
plot(time,signal,'k')
xlabel('Time'), ylabel('Amplitude')

%% "Static" Power Spectrum from FFT
% Vector of frequencies
hz = linspace(0,srate/2,floor(npnts/2)+1);

% Compute signal power
sigpower = abs(fft(signal)/npnts).^2;

% Plot the power spectrum
figure(2)
subplot(211)
plot(hz,sigpower(1:length(hz)),'k','linew',2)
set(gca,'xlim',[0 max(frex)+10])
xlabel('Frequency (Hz)'), ylabel('Power')
title('Power spectrum from FFT')

% Show dashed vertical lines with simulated frequencies
hold on
plot(repmat(frex,2,1),repmat(get(gca,'ylim')',1,length(frex)),'r--');

%% Welch's Method
% Window length in time points
winlength = 1000;

% Window onset times
winonsets = round(linspace(1,npnts-winlength,50));

% Different-length signal needs a different-length Hz vector
hzW = linspace(0,srate/2,floor(winlength/2)+1);

% Initialize the power matrix 
sigpower2 = zeros(length(winonsets),length(hzW));

% Loop over frequencies
for wi=1:length(winonsets)
    % Get a chunk of data from this time window
    datachunk = signal(winonsets(wi):winonsets(wi)+winlength);
    
    % Compute its power
    tmppow = abs(fft(datachunk)/winlength).^2;
    
    % Enter into matrix
    sigpower2(wi,:) = tmppow(1:length(hzW));
end

% Plot the power spectrum
subplot(212)
plot(hzW,mean(sigpower2,1),'k','linew',2)
set(gca,'xlim',[0 max(frex)*2])
xlabel('Frequency (Hz)'), ylabel('Power')
title('Power spectrum from Welch''s method')

% Show  simulated frequencies
hold on
plot(repmat(frex,2,1),repmat(get(gca,'ylim')',1,length(frex)),'r--');

%% end