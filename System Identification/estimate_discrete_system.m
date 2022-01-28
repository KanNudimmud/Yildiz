%% Estimation of a Discrete-Time System with 2 Delays,Zeros & 3 Poles
%% Frequency Response
% Create a continuous-time system
Npts = 512; T=.5; t=(0:1:Npts-1)*T;
RandStream.setGlobalStream(RandStream('mcg16807','Seed',15));

nn = [1 6 -1/6]; dd=[1 1 8 -0.1];
H = tf(nn,dd,'InputDelay',0.7,'OutputDelay',1.7)

% Convert to discrete-time
[num,den] = c2dm(H.Numerator,H.Denominator,T,'zoh');

w_arx=logspace(-2,0,300)*pi/T;

% Compute input sequence for identification
u = idinput(Npts,'prbs'); % input signal for identification
yc = dlsim(num,den,u); % compute "clean" output of G
v = randn(Npts,1); % noise - note that this is white, gaussian
LL = 0.25*(yc'*yc)/(v'*v); % scale so energy in sensor noise 1/4 times
v = sqrt(LL)*v; % the energy in the "clean" signal
y = yc+v; % actual output y=Gu+v
Z = iddata(y,u,T); % data available to identification
% Compute input sequence for validation
u_val = idinput(Npts,'prbs'); % input signal for validation
yc_val = dlsim(num,den,u_val); % "clean" output of G
v_val = sqrt(LL)*randn(Npts,1); % noise - white, gaussian
y_val = yc_val+v_val; % actual output y=Gu+v
Z_val = iddata(y_val,u_val,T); % data available for validation
% Plot portions of input and output signals after initial transient
figure; plot(t,[u v]); legend('u[k]','v[k]');
title('Inputs to the system being identified');
xlabel('Time (s)'); ylabel('Amplitude (unitless)');

figure; plot(t,[yc y]); legend('"clean" y_c[k]','measured y[k]');
title('Outputs from the system being identified');
xlabel('Time (s)'); ylabel('Amplitude (unitless)');

%% Identify System Models 
% Frequency response of actual system, and "SPA" frequency resp model
[mag,ph,w]=dbode(num,den,T,w_arx); % get "true" magnitude and phase resp
G = spa(Z,64,w,[],T); [amp,phas,w]=bode(G); w = squeeze(w);
amp = squeeze(amp); phas = squeeze(phas);

% ARX model with na=2; nb=2; nk=1 (ARX221)
M_arx221 = arx(Z,'na',2,'nb',2,'nk',1);
[m_arx221,p_arx221,w_arx221]=bode(M_arx221); w_arx221 = squeeze(w_arx221);
m_arx221 = squeeze(m_arx221); p_arx221 = squeeze(p_arx221);
[a_arx221,b_arx221,c_arx221,d_arx221,f_arx221] = polydata(M_arx221);

% ARX model with na=4; nb=4; nk=1 (ARX441)
M_arx441 = arx(Z,'na',4,'nb',4,'nk',1);
[m_arx441,p_arx441,w_arx441]=bode(M_arx441); w_arx441 = squeeze(w_arx441);
m_arx441 = squeeze(m_arx441); p_arx441 = squeeze(p_arx441);
[a_arx441,b_arx441,c_arx441,d_arx441,f_arx441] = polydata(M_arx441);

% ARMAX model with na=2; nb=2; nc=2; nk=1 (ARMAX2221)
M_armax=armax(Z,'na',2,'nb',2,'nc',2,'nk',1);
[a_armax,b_armax,c_armax,d_armax,f_armax]=polydata(M_armax);
[m_armax,p_armax,w_armax]=bode(M_armax); w_armax = squeeze(w_armax);
m_armax = squeeze(m_armax); p_armax = squeeze(p_armax);

% Box-Jenkins model with nb=2; nc=2; nd=2; nf=2; nk=1 (BJ22221)
% y(t) = [B(q)/F(q)] u(t-nk) + [C(q)/D(q)] e(t)
M_bj=bj(Z,'nb',2,'nc',2,'nd',2,'nf',2,'nk',1);
[m_bj,p_bj,w_bj]=bode(M_bj); w_bj = squeeze(w_bj);
m_bj = squeeze(m_bj); p_bj = squeeze(p_bj);
[a_bj,b_bj,c_bj,d_bj,f_bj]=polydata(M_bj);

% OE model with nb=2; nf=2; nk=1;
M_oe = oe(Z,'nb',2,'nf',2,'nk',1);
[m_oe,p_oe,w_oe]=bode(M_oe); w_oe = squeeze(w_oe);
m_oe = squeeze(m_oe); p_oe = squeeze(p_oe);
[a_oe,b_oe,c_oe,d_oe,f_oe]=polydata(M_oe);

% Now, plot Bode plots
figure; loglog(w,mag,w,amp,w_arx221,m_arx221,w_arx441,m_arx441);
title('Bode mag. plots of several system id models');
ylabel('Magnitude'); xlabel('Frequency (rad s^{-1})');
legend('Actual','SPA','ARX221','ARX441'); axis([.08 8 1e-2 5]);

figure; semilogx(w,ph,w,phas,w_arx221,p_arx221,w_arx441,p_arx441);
title('Bode phase plots of several system id models');
xlabel('Frequency (rad s^{-1})'); ylabel('Phase (deg)');
legend('Actual','SPA','ARX221','ARX441'); axis([.08 8 -270 0]);

figure; loglog(w,mag,w_oe,m_oe,w_armax,m_armax,w_bj,m_bj);
title('Bode mag. plots of several system id models');
ylabel('Magnitude'); xlabel('Frequency (rad s^{-1})');
legend('Actual','OE221','ARMAX2221','BJ22221'); axis([.08 8 1e-2 5]);

figure; semilogx(w,ph,w_oe,p_oe,w_armax,p_armax,w_bj,p_bj);
title('Bode phase plots of several system id models');
xlabel('Frequency (rad s^{-1})'); ylabel('Phase (deg)');
legend('Actual','OE221','ARMAX2221','BJ22221'); axis([.08 8 -270 0]);

%% Unit-Pulse Response
% True and estimated system discrete-time unit-pulse responses
Ntime=30;
y_act = dimpulse(num,den,Ntime);
y_arx221 = dimpulse(b_arx221,a_arx221,Ntime);
y_arx441 = dimpulse(b_arx441,a_arx441,Ntime);
y_armax = dimpulse(b_armax,a_armax,Ntime);
y_oe = dimpulse(b_oe,f_oe,Ntime);
y_bj = dimpulse(b_bj,f_bj,Ntime);

figure; stem([0:Ntime-1]*T,[y_act y_arx221 y_arx441],'filled'); hold on
stem([0:Ntime-1]*T,y_act,'filled');
legend('True system','ARX221','ARX441');
title('Discrete impulse responses')
xlabel('Time (sec)'); ylabel('Output amplitude')

figure; stem([0:Ntime-1]*T,[y_act y_armax y_oe y_bj],'filled'); hold on
stem([0:Ntime-1]*T,y_act,'filled');
legend('System','ARMAX2221','OE221','BJ22221');
title('Discrete impulse responses')
xlabel('Time (sec)'); ylabel('Output amplitude')

% Model validation
% Compute residuals
e_arx2 = resid(M_arx221,Z); e_arx2 = e_arx2.OutputData;
e_arx4 = resid(M_arx441,Z); e_arx4 = e_arx4.OutputData;
e_arm = resid(M_armax,Z); e_arm = e_arm.OutputData;
e_bj = resid(M_bj,Z); e_bj = e_bj.OutputData;
e_oe = resid(M_oe,Z); e_oe = e_oe.OutputData;
mean([v e_arx2 e_arx4 e_arm e_bj e_oe])'

% Plot a histogram of the residuals for ARX221
figure; hist([e_arx2 v],-2:0.2:2); axis([-1.6 1.6 0 160])
title('Residual histogram for ARX221');
ylabel('Count');xlabel('Value of residual')
legend('Model fit','Actual');

% Similar for the other cases. Omitting code for figure formatting...
figure; hist([e_arx4 v],-2:0.2:2);axis([-1.6 1.6 0 160])
title('Residual histogram for ARX441');
figure; hist([e_arm v],-2:0.2:2);axis([-1.6 1.6 0 160])
title('Residual histogram for ARMAX2221');
figure; hist([e_bj v],-2:0.2:2);axis([-1.6 1.6 0 160])
title('Residual histogram for BJ22221');
figure; hist([e_oe v],-2:0.2:2);axis([-1.6 1.6 0 160])
title('Residual histogram for OE221');

%% Model Validation using Correlations
% Create new figures; call "resid" with no outputs to plot residuals
figure; resid(M_arx221,Z);
figure; resid(M_arx441,Z); figure; resid(M_armax,Z);
figure; resid(M_bj,Z); figure; resid(M_oe,Z);

%% end.