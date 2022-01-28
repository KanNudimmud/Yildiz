%% Hydraulic System
% Set Parameters
del_p  = 70;
r_s    = (del_p/5)*60e-3;
l_p    = 20e-3;
r_p    = 0.203;
rho_g  = 1e4;
D      = 1e-2;
A      = pi*D^2/4;
e_h    = (rho_g/A) * (1e-6/132);
w_n_sq = 2*e_h/l_p;
tau_p  = r_p/l_p;

%%  end.