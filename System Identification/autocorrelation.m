%% Autocorrelation of White Noise
% System Identification toolbox is needed.
clear all
clc
%% Create a system
% System parameters
var   = 0.5; % variance of the input noise
n     = 1024; % number of data points
alpha = 2; % input scaling

% Unit impulse input
uPulse = alpha * eye(n,1); % eye fnc. creates identity matrix

% Polynomial parameters
x1 = -0.8; 
y1 = 1;
z1 = 1;

% Transfer function
tf1 = idpoly([1 x1], [0 0 y1], [z1], [1]); % creates polynomial with parameters

% Error
err = sqrt(var) * randn(n,1);

% Output without error
y = sim(tf1, uPulse, 0*err); % simulate dynamic system

% Input and output with error (white noise)
uWhite = alpha * 5* randn(n,1);
yWhite = sim(tf1, uWhite, err); 

%% Calculate Correlation
% Correlation between input and output
[Ryu, lagsyu]= xcorr (yWhite, uWhite,'biased');

% Find which correlation measured at last
indyu=find(lagsyu==0);

% Correlation for input
[Ruu, lagsuu]=xcorr(uWhite, uWhite,'biased');

% Find which correlation measured at last
induu=find(lagsuu==0);

% Create symmetric Toeplitz matrix (used to model LTI systems)
RuMat=toeplitz(Ruu(induu:end)); 

% Pulse response
ghat=RuMat \ Ryu(induu:end);

%% Display Results
% Initialize x-axis
k=lagsyu(indyu:end);

% Display pulse response
figure (1)
stem (0:19, y(1:20)/alpha)
hold on;
stem (k, ghat,'r>')
legend ('Actual','Correlation approx.')
title ('Pulse response estimate')
ylabel ('Pulse resp est:\itu\rm[\itk\rm] white, correlation method')
xlabel ('Time (Samples')
xlim([0 20])

% Display autocorrelation
figure (2)
stem (k, Ruu(induu:end),'rs')
title ('Approx. autocorrelation of white-noise input')
ylabel ('Approximate autocorrelation')
xlabel('Delay (Samples)')
xlim([0 20])

%% end