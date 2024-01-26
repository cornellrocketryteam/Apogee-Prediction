% apogee prediction using Kalman filter

% xdot = Ax + bu
% xdot_est = Ax_est + Bu + L(y-Cx_est)
% y = Cx
% u = -Kx + k_r*r