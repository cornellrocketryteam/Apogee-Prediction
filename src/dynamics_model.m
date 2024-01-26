clc
clear
close all

data = readtable('Flight_Data.csv');

time = table2array(data(:,1));
altitude = table2array(data(:,2));
velocity = table2array(data(:,4));
flight = table2cell(data(:,9));

L = length(time);

time1 = [];
time2 = [];
time3 = [];
time4 = [];

w = 1;
x = 1;
y = 1;
z = 1;

for i = 1:L
    launch = flight(i);
    if char(launch) == 'TL 22'
        time1(w) = time(i);
        altitude1(w) = altitude(i);
        rrc3_v1(w) = velocity(i);
        w = w + 1;
    elseif char(launch) == 'CL 22'
        time2(x) = time(i);
        altitude2(x) = altitude(i);
        rrc3_v2(x) = velocity(i);
        x = x + 1;
    elseif char(launch) == 'TL 23'
        time3(y) = time(i);
        altitude3(y) = altitude(i);
        rrc3_v3(y) = velocity(i);
        y = y + 1;
    else
        time4(z) = time(i);
        altitude4(z) = altitude(i);
        rrc3_v4(z) = velocity(i);
        z = z + 1;
    end

end

R = 287; % ideal gas constant of air
g0 = 9.8; % gravitational constant
a = -0.0065; 
rho_ref = 1.225; 
T_ref = 288.15;

% densities for altitudes
rho1 = (1/16.0185)*(1+(a*(altitude1*0.3048)/T_ref)).^((-g0/(a*R))-1); % ft/lb^3
rho2 = (1/16.0185)*(1+(a*(altitude2*0.3048)/T_ref)).^((-g0/(a*R))-1);
rho3 = (1/16.0185)*(1+(a*(altitude3*0.3048)/T_ref)).^((-g0/(a*R))-1);
rho4 = (1/16.0185)*(1+(a*(altitude4*0.3048)/T_ref)).^((-g0/(a*R))-1);

% distance to apogee
d1 = 10700 - altitude1; % ft
d2 = 10345 - altitude2;
d3 = 12315 - altitude3;
d4 = 10067 - altitude4;

TL_t = [time1 time3]; % single matrix for test launch times
CL_t = [time2 time4]; % single matrix for comp launch times
TL_d = [d1 d3]; % single matrix for test launch distances to apogee
CL_d = [d2 d4]; % single matrix for comp launch distances to apogee
TL_rrc3v = [rrc3_v1 rrc3_v3]; % filtered velocity from test launch rrc3's
CL_rrc3v = [rrc3_v2 rrc3_v4]; % filtered velocity from comp launch rrc3's

% TIME FOR SOME DYNAMICS
t = linspace(8.05,27,1000);
URRG = 728; % URRG start in ft
spaceport = 4595; % spaceport start in ft

[t,x1] = ode45(@(t,X) pic(t,X,URRG),t,[altitude1(1);rrc3_v1(1)]);
[t,x2] = ode45(@(t,X) pic(t,X,spaceport),t,[altitude2(1);rrc3_v2(1)]);




figure 

scatter(time1,altitude1,'g')
hold on
plot(t,x1(:,1),'k')

scatter(time2,altitude2,'r')
plot(t,x2(:,1),'k')
xlabel('Time (s)','FontSize',16)
ylabel('Altitude (ft)','FontSize',16)
title('Post-Burnout RRC3 Altitude Data for CRT Flights','FontSize',20)
hold off


function Xdot = pic(t,X,xstart)
    x = X(1);
    xdot = X(2);
    R = 287;
    g0 = 9.81;
    g = 32.17; % ft/s
    a = -0.0065;
    rho_ref = 1.225;
    T_ref = 288.15;
    Cd = 0.536; % coefficient of drag
    S = pi*0.25^2; % cross sectional area of rocket in ft^2
    m = 140;

  

    beta = 0.5*Cd*S;
   

    rho = (1/16.0185)*(rho_ref*(1+(a*((x+xstart)*0.3048)/T_ref)).^((-g0/(a*R))-1)); % ft/lb^3

    xdoubledot = -g - beta/m*rho*xdot^2;
    Xdot = [xdot; xdoubledot];
    
end


