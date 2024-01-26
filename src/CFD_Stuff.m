clc
clear
close all

mach = [0.6045,0.5614,0.5182,0.4750,0.4318,0.38886,0.3455,0.3023]';
F_AB = 2*[145.17,125.79,106.3,88.66,72.45,58.65,45.91,34.81]';
F_LV = 2*[95.7,83.97,72.02,60.98,50.83,41.92,33.42,26.02]';

F_total = F_AB + F_LV;

p = polyfit(mach,F_total,2);
m = linspace(0.1,0.7,1000);

theoretical = p(1)*m.^2 + p(2)*m.^1 + p(3);


scatter(mach,F_total,'k')
hold on
plot(m,theoretical,'r','LineWidth',1.5)
xlabel('Mach Number','FontSize',16)
ylabel('Drag Force on LV (N)','FontSize',16)
title('Drag at Varying Mach Numbers with Airbrakes Deployed','FontSize',16)

R = 287;
g0 = 9.81;
g = 32.17; % ft/s
a = -0.0065;
rho_ref = 1.225;
T_ref = 288.15;

xstart = 4595;

% 6426 ft at start of comp launch 2022 data
% 510.29 ft/s at start of comp launch 2022 data
% 10310 apogee in real life

