clc
clear
close all

% Test launch predictions are too high. This indicates we are
% underestimating the drag force at low altitudes. Air density is too small

% Comp launch predictions are too low. This indicates we are overestimating
% the drag force at high altitudes. Air density is too big.

% Repair model so density needs to be bigger at low altidues and smaller at high
% altitudes



R = 287;
g0 = 9.81;
g = 32.17; % ft/s
a = -0.0065;
rho_ref = 1.225;
T_ref = 350;
x = linspace(0,15000,15001);
rho = (1/16.0185)*(rho_ref*(1+(a*((x)*0.3048)/T_ref)).^((-g0/(a*R))-1)); % ft/lb^3


plot(x,rho,'b')
xlabel('Altitude above sea level [lb]')
ylabel('Air Density [lb/ft^3]')
ylim([0,0.15])
title('Air Density by Altitude above Sea Level','FontSize',16)


hold on

R = 85;
g0 = 9.81;
g = 32.17; % ft/s
a = -0.0065;
rho_ref = 2.5;
T_ref = 300;
rho2 = (1/16.0185)*(rho_ref*(1+(a*((x)*0.3048)/T_ref)).^((-g0/(a*R))-1)); % ft/lb^3

plot(x,rho2,'r')

R = 287;
g0 = 9.81;
g = 32.17; % ft/s
a = -0.0065;
rho_ref = 1.225;
T_ref = 350;
x = linspace(0,15000,15001);
rho3 = (1/16.0185)*(rho_ref*(1+(a*((x)*0.3048)/T_ref)).^((-g0/(a*R))-1)); % ft/lb^3

plot(x,rho3,'g')

legend('Original Density Model','Adjusted Density Model','Adjusted Density due to T')
hold off
