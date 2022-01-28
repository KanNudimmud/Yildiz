%% Dryer Example (System Identification Toolbox needed)
% Load the data - input is current, output is temperature-
load dryer2;

% Scale the data
u2 = u2 - mean(u2);
y2 = y2 - mean(y2);

% Parameters
Ts = 0.08;
x = -length(u2)+1:length(u2)-1;

% Cross-Correlation for delay
Ryu = xcorr(y2,u2);

figure(1)
stem(x,Ryu,'filled')
title('Cross-correlation between y and u')
xlabel('Time shift(samples)'), ylabel('Correlation')

figure(2)
stem(x,Ryu,'filled')
title('Cross-correlation between y and u')
xlabel('Time shift(samples)'), ylabel('Correlation')
xlim([0,10])

% Split the data to train and test
ue = u2(1:length(u2)/2);
ye = y2(1:length(y2)/2);

uv = u2(length(u2)/2+1:end);
yv = y2(length(u2)/2+1:end);

% Encapsulate splitted input and output respect to sampling rate
ze = iddata(ye,ue,Ts);
zv = iddata(yv,uv,Ts);

% Set names of parameters
ze.InputName = 'Current';
zv.InputName = 'Current';
ze.InputUnit = 'A';
zv.InputUnit = 'A';

ze.OutputName = 'Temperature';
zv.OutputName = 'Temperature';
ze.OutputUnit = '\circC';
zv.OutputUnit = '\circC';

ze.TimeUnit = 'sec';
zv.TimeUnit = 'sec';

% Plot the data
figure(3)
plot(ze)

% Step Response
figure(4)
step(impulseest(ze),2)
xlim([-1 2])

% Estimate time delay
nk = delayest(ze);

% Bode Diagram
Ge = spa(ze);

figure(5)
bode(Ge)

% Estimate Frequency Response
z       = iddata(y2,u2,Ts);
[m,p,w] = bode(spa(z));
w       = squeeze(w);
m       = squeeze(m);

figure(6)
semilogx(w,20*log10(m))
title('Estimated Frequecny Response')
xlabel('Frequency(Hz)')
ylabel('Magnitude(dB)')

% Generate a matrix for ARX model (combinations)
NN = struc(1:5,1:5,nk);

% Loss function as a vector
V = arxstruc(ze,zv,NN);

% Select best combination
theStruc = selstruc(V);

% Estimate parameters of ARX model
Marx = arx(ze,theStruc);
present(Marx)

% Find polynomial roots
marxPoles = roots(Marx.a);

% Compare the original data and estimated
figure(7)
compare(zv,Marx)

% Compute residuals
figure(8)
resid(zv,Marx)

% Compare original and predicted data
yp = predict(Marx,zv,10);

figure(9)
plot(Ts*(1:length(ue)),[zv.OutputData,yp.OutputData])
title('Comparing system output to predicted output')
xlabel('Time(seconds)')
ylabel('Temperature(C)')
legend({'Actual data';'Predicted data'})

% Compute the prediction error
err = pe(Marx,zv);

figure(10)
h   = bodeplot(spa(err,[],logspace(-2,2,200)));
showConfidence(h,3)
title('Power Spectrum for Temperature')

%% end.