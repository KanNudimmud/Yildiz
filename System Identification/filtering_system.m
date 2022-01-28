%% Mass-Spring-Damper System   
clear all,clc
%% Initial Values
x_0=0;v_0=0;
syms s;
m=2;c=0.4;k=0.6;
F_s=10*1/s;
X_s=F_s*1/(m*s^2+c*s+k);
X_t_cont=ilaplace(X_s);

% Display
figure(1),fplot([X_t_cont])
xlim([0 60]),ylim([0 30])
grid on
xlabel('Time (Sec)'),ylabel('Position (m)')

%% Initial Parameters
dt = 0.001;ti=0;tf=15;
t=[ti:dt:tf]';
length_of_loop=(tf-ti)/dt;

% Forward difference
X_t_f=zeros(length_of_loop,1);
V_t_f=zeros(length_of_loop,1);

% input
F = cos(2*t);

clear s;
s.A=[1,dt;
    -k*dt/m,(1-dt*c/m)];

% Measurement noise 
MNstd=0.004;MNV=MNstd*MNstd;

% Process noise
PNstd=0.002;PNV=PNstd*PNstd;

s.Q=1e-12*eye(2)*PNV;

% Define measurement function to  return the state
s.H=[1,0;
    0,1];
C_out=[1,0;
    0,1];

% Define a measurement error
s.R=0.00001*eye(2)*MNV

% Use control to include gravity
s.B=[0;
    dt/m];

% Initial state
s.x=[x_0,v_0]';
s.P=eye(2)*MNV;
s.detP=det(s.P) % keep track of the noise
s.z=zeros(2,1);

% Create dynamics matrices
tru = zeros(length_of_loop,2) % kalman dynmacis
trutr =zeros(length_of_loop,2); % norm dynamics
tru_noisy =zeros(length_of_loop,2); % noisy dynamic
tru_noisy_out =zeros(length_of_loop,2); % noisy out dynamic
tru_output = zeros(length_of_loop,2);

tru(1,:)=[x_0 v_0];
detP(1,:)=s.detP;

Atr = [1, dt;
      -k*dt/m, (1 - dt*c/m)];
Btr = [0;
      dt/m];

for i = 1:1:length_of_loop-1
    noise_process = PNstd*randn(2,1);
    noise_measurement = MNstd*randn(2,1);
    % true system dynamic
    trutr(i+1,:) = Atr*trutr(i,:)' + Btr*F(i);
    % noisy dynamic
    tru_noisy(i+1,:) = Atr*tru_noisy(i,:)' + Btr*F(i) + noise_process;
    tru_noisy_out(i+1,:) = C_out*tru_noisy(i+1,:)' + noise_measurement;
    % kalman dynamic
    tru(i+1,:)=s(i).A*tru(i,:)'+ s(i).B*F(i) + noise_process;
    s(i).z = s(i).H*tru(i+1,:)' + noise_measurement;
    s(i+1) = kalmanf(s(i),F(i));
    detP(i+1) = s(i+1).detP;
    % true output of kalman filter
    tru_output(i+1,:) = s(i+1).x;
end

figure(2)
plot(t(1:end-1),trutr(:,1))
hold on
plot(t(1:end-1),tru_noisy_out(:,1))
hold on
plot(t(1:end-1),tru_output(:,1))
legend('True','Noise','Kalman')

Meas_err_noise = trutr(:,1) - tru_noisy_out(:,1);
Meas_Err_noise_cov = sum(Meas_err_noise.*Meas_err_noise)/length(Meas_err_noise);
Meas_err_kalman = trutr(:,1) - tru_output(:,1);
Meas_Err_kalman_cov = sum(Meas_err_kalman.*Meas_err_kalman)/length(Meas_err_kalman);

%% end.