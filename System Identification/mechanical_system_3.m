%% Create a dynamic system (mechanical system with 3 mass,damper,spring)
% Initialize parameters 
m1 = 4; m2 = 3; 
c1 = 2; c2 = .4; c3 = .7; 
k1 = 20; k2 = 30; k3 = 40;

% Set-up matrices
M = diag([m1 m2 2*m2]);
K = [k1+k2 -k2 0; -k2 k2+k3 -k3; 0 -k3 k3];
C = [c1+c2 -c2 0; -c2 c2+c3 -c3; 0 -c3 c3];
a = [zeros(3) eye(3); -M\[K C]]; % state matrix
b = [zeros(3,1); M\[1;0;0]]; % input matrix
c = [m1 0 m2 0 0 0]; % output matrix
d = 0; % direct transmission matrix

% Sampling parameters
Npnts = 1024; Ts = .2; Tf = (Npnts-1)*Ts; t = (0:Ts:Tf)';

%% Continuous Time-domain Response of the System to Unit Impulse Input
sys = ss(a,b,c,d);
figure(1)
impulse(sys)
xlim([0 120])

%% Discrete-time Transfer Function of the System
% Convert continuous-time system to discrete-time 
sysd = c2d(sys,Ts);

% Simulate the system
figure(2)
lsim(sysd,ones(1024,1),t)
xlim([0 120])

%%  Discrete Time-domain Response of the System to Unit Impulse Input
figure(3)
impulse(sysd)
xlim([0 120])

%% Discrete Time-domain Response of the System to Random Input in the presence of White Sensor Noise
% Create noise
sigma = 1/8; % noise size
u     = idinput(length(t),'prbs'); % pseudo random noise input
w     = [u sigma*randn(size(t))]; % sensor noise
w     = w(:,1); % take one channel for consistency

% Simulate the system
figure(4)
lsim(sysd,w,t)
xlim([0 120])

%% Estimate the system 
% ETFE
y=lsim(sys,w,t);
z=[y u];

% Filtering of both y and u to chop out HF stuff
zf=idfilt(z,4,.8);

f  = logspace(-2,1,400); 
g  = freqresp(a,b,c,d,1,1j*2*pi*f);
phg= unwrap(angle(g))*180/pi;

% ETFE for different M_inputs
Minps = [512,256,128,64];
for theMinp = 1:length(Minps)
    Minp = Minps(theMinp);
    ghat = etfe([y u],Minp,Npnts/2,Ts);
    [what,ghm,ghp]=getff(ghat,1,1);
    figure(4+theMinp)
    subplot(121)
    loglog(2*pi*f,abs(g),what,ghm,'.')
    xlabel('Frequency (rad s^{-1})'), ylabel('Magnitude')
    title(sprintf('Mag. est.: ETFE with N=%d, M_{input} = %d',Npnts,Minp))
    xlim([.1 20]),ylim([.001 10])
    subplot(122)
    semilogx(2*pi*f,phg,what,ghp,'.')
    xlabel('Frequency (rad s^{-1})'), ylabel('Phase (deg)')
    title(sprintf('Phase est: ETFE with N=%d, M_{input} = %d',Npnts,Minp))
    xlim([.1 20]),ylim([-500 100])
end

% SPA for different gammas
gammas = [128,256,512];
for theGamma = 1:length(gammas)
    winGam=gammas(theGamma);
    ghats=spa([y u],winGam,what,[],Ts);
    [what,ghms,ghps]=getff(ghats,1,1);
    figure(8+theGamma)
    subplot(121)
    loglog(2*pi*f,abs(g),what,ghms,'.')
    xlabel('Frequency (rad s^{-1})'), ylabel('Magnitude')
    title(sprintf('Mag. est.: SPA with N=%d, \\gamma = %d',Npnts,winGam))
    xlim([.1 20]),ylim([.001 10])
    subplot(122)
    semilogx(2*pi*f,phg,what,ghps,'.')
    xlabel('Frequency (rad s^{-1})'), ylabel('Phase (deg)')
    title(sprintf('Phase est.: SPA with N=%d, \\gamma = %d',Npnts,winGam))
    xlim([.1 20]),ylim([-500 100])
end

%% end