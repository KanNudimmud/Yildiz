%% Convolution of two signals
clear all, clc

% Create constants
a = linspace(-10,10,100);
b = linspace(-1,1,100);
t = linspace(0,3,100);
A = linspace(-6,6,100);
f = 512;
w = 2*pi*f;

% Initialize signals
f =  a .* exp(-b.*t);
g =  A .* sin(w*t);

%% Convolution with a built-in function
C = conv2(f,g,'same');

% Display functions and their convolution
figure(1),clf
plot(t,f, t,g, t,C)
legend({'a*e^(-bt)';'A*sin(wt)';'Conv2'})

%% Convolution in Time Domain 
% Convolution sizes
n = length(f); % f and g have same sizes
nConv = 2*n - 1;

half = floor(n/2);

% flipped version of g
gflip = g(end:-1:1);

% zero-padded data for convolution
data = [ zeros(1,half) f zeros(1,half) ];

% initialize convolution output
Ctime = zeros(1,nConv);

% run convolution
for i=half+1:nConv-half
    % get a chunk of data
    chunk = data(i-half+1:i+half);
    
    % compute dot product
    Ctime(i) = sum( chunk.*gflip );
end

% cut off edges
Ctime = Ctime(half+1:end-half);

%% Convolution in Frequency Domain 
% spectra of f and g
fSpec = fft(f,nConv);
gSpec = fft(g,nConv);

% element-wise multiply
fXg = fSpec .* gSpec;

% inverse FFT to get back to the time domain
Cfreq = ifft( fXg );

% cut off edges
Cfreq = Cfreq(half+1:end-half);

%% Visualize the Comparison
figure(2)
subplot(311)
plot(C,'r-'),legend('Built-in')
subplot(312)
plot(Ctime,'b-'),legend('Time domain')
subplot(313)
plot(Cfreq,'k-')
legend('Freq. domain')

%% end