%% Create a dynamic system
% Initialize parameters 
m1 = 4; m2 = 3; c1 = 2; c2 = .4; c3 = .7; k1 = 20; k2 = 30; k3 = 40;

% Set-up matrices
M = diag([m1 m2 2*m2]);
K = [k1+k2 -k2 0; -k2 k2+k3 -k3; 0 -k3 k3];
C = [c1+c2 -c2 0; -c2 c2+c3 -c3; 0 -c3 c3];
a = [zeros(3) eye(3); -M\[K C]]; % state matrix
b = [zeros(3,1); M\[1;0;0]]; % input matrix
c = [m1 0 m2 0 2*m2 0]; % output matrix
d = 0; % direct transmission matrix

% Sampling parameters
Npnts = 1024; Ts = .2; Tf = (Npnts-1)*Ts; t = (0:Ts:Tf)';

% Create noise
sigma = 1/8; % noise size
u     = idinput(length(t),'prbs'); % pseudo random noise input
w     = [u sigma*randn(size(t))]; % sensor noise

% Simulate the system
y  = lsim(a,[b 0*b],c,[d 1],w,t);
z  = [y u];
zf = idfilt(z,4,.8);

%% Display Frequency Response
% Initialize parameters
f   = logspace(-2,1,400); 
g   = freqresp(a,b,c,d,1,1j*2*pi*f);
phg = unwrap(angle(g)) * 180/pi;

% Visualize true frequency response
figure(1), loglog(2*pi*f,abs(g));
xlabel('Frequency (rad s^{-1})'), ylabel('Magnitude')
title('True Frequency Response of the System')
xlim([.1 50]),ylim([.001 5])

% Visualize sample time sequence
figure(2)
plot(t(1:101),z(1:101,1),t(1:101),w(1:101,2),'--')
xlabel('Time (s)'), ylabel('Output')
title('Sample Time Sequence: Indicates SNR')
legend(['y';'v']),ylim([-1 1])

%% ETFE
Minps = [512,256,128,64];
for Mi = 1:length(Minps)
    Minp=Minps(Mi);
    ghat=etfe([y u],Minp,Npnts/2,Ts); 
    [what,ghm,ghp]=getff(ghat,1,1);
    figure(2+Mi),subplot(121)
    loglog(2*pi*f,abs(g),what,ghm,'.');
    xlabel('Frequency (rad s^{-1})'), ylabel('Magnitude')
    title(sprintf('Mag. est.: ETFE with N=%d, M_{input} = %d',Npnts,Minp))
    subplot(122)
    semilogx(2*pi*f,phg,what,ghp,'.');
    xlabel('Frequency (rad s^{-1})'),ylabel('Phase (deg)')
    title(sprintf('Phase est: ETFE with N=%d, M_{input} = %d',Npnts,Minp))
end

%% SPA
gammas = [128,256,512];
for Gi = 1:length(gammas)
    winGam=gammas(Gi);
    ghats=spa([y u],winGam,what,[],Ts); 
    [what,ghms,ghps]=getff(ghats,1,1);
    figure(6+Gi),subplot(121)
    loglog(2*pi*f,abs(g),what,ghms,'.');
    xlabel('Frequency (rad s^{-1})'), ylabel('Magnitude')
    title(sprintf('Mag. est.: SPA with N=%d, \\gamma = %d',Npnts,winGam))
    subplot(122)
    semilogx(2*pi*f,phg,what,ghps,'.');
    xlabel('Frequency (rad s^{-1})'), ylabel('Phase (deg)')
    title(sprintf('Phase est.: SPA with N=%d, \\gamma = %d',Npnts,winGam))
end

%% end.