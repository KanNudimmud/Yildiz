%% Kalman Filtering on Mechanical System with 2 Mass-Damper & 3 Spring
%% Model Construction
% Initialize parameters 
m1 = 1; m2 = 2; 
c1 = 0.01; c2 = .05; 
k1 = 2; k2 = 3; k3 = 1;

M = diag([m1,m2]);
k = [k1+k2 -k2; -k2 k2+k3];
c = [c1+c2 -c2; -c2 c2];
A = [zeros(2) eye(2); -M\[k c]];
B = [zeros(2); inv(M)];
Cme = [0 1 0 0]; % assume only the displacement of the second mass object can be measured
C = eye(4); % to extract all the states for result comparison
D = [0 0];
Csys = ss(A,B,Cme,0); % original system (continuous system)

% Discretise the continuous system 'sys_c'
Ts = 0.1; Dsys = c2d(Csys,Ts);

% C matrix remains the same
A = Dsys.A; B = Dsys.B;

% External forces
t = [0:Ts:60];
omega_f = 1; mag = 2;
f1 = mag*sin(omega_f*t);
f2 = zeros(1,length(t));
F = [f1; f2];

% Initial condition
x0 = [0;0;0;0];

%% Information for the Noises
% The variances for the process noise and the measurement noise
covq = [1]; covr = [0.01];

% Process noise array
w = sqrt(covq)*randn(2,length(t));

% Measurement noise array
v = sqrt(covr)*randn(1,length(t));

%% Get the Measurement 
% Simulation for measurement 'y_o', which is contaminated by both process
% noise and measurement noise. This measurement 'y_o' will be provided to
% the Kalman filter
sys_n = ss(A,[B B zeros(4,1)],Cme,[zeros(1,2) zeros(1,2) 1],Ts);
y_o = lsim(sys_n,[F;w;v],t,x0);

%% Simulate True State Vector
% Simulation for the true states (three of them unmeasureable),
% this system contains the process noise but no measurement noise
sys_r = ss(A,[B B zeros(4,1)],C,0,Ts);
y_r = lsim(sys_r,[F;w;v],t,x0);

%% Design Kalman Filter
y_e = zeros(length(t),1);
L = zeros(4,length(t));
x = x0;
P = B*(covq.*eye(2))*B';
for k_time = 1:length(t)
  % Measurement update
  Mxn = P*Cme'/(Cme*P*Cme'+covr);
  x = x + Mxn*(y_o(k_time)-Cme*x);   % x[n|n]
  P = (eye(4)-Mxn*Cme)*P;     % P[n|n]
  y_e(k_time) = Cme*x;
  L(:,k_time) = Mxn;

  % Time update
  x = A*x + B*F(:,k_time);        % x[n+1|n]
  P = A*P*A' + B*(covq.*eye(2))*B';     % P[n+1|n] 
end

%% Plot the Results
figure(1)

% Plot the reference, observation and estimation
subplot(211)
plot(t,y_r(:,2),'LineWidth',1,'Color','r'),hold on
plot(t,y_o,'Marker','*','Color','b','LineStyle','none')
plot(t,y_e,'LineWidth',1,'Color','g'), hold off
%axis([0 20 -2 2])
xlabel('Time (seconds)')
ylabel('Displacement (meters)')
title('Displacement for the second mass')
legend('Reference', 'Measurement', 'Estimation')

% Plot the error
subplot(212)
plot(t,y_o-y_r(:,2),'LineWidth',1,'Color','b'), hold on
plot(t,y_e-y_r(:,2),'LineWidth',1,'Color','g')
xlabel('Time (Seconds)')
ylabel('Error')
title(['Measurement error=',num2str(sum(abs(y_o-y_r(:,2)))/length(y_o)),'    Estimation error=',num2str(sum(abs(y_e-y_r(:,2)))/length(y_o))])
legend('Measurement', 'Estimation')

%% end