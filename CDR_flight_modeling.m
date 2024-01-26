clc
close all

velocity_smoothing

% elimate a bunch of zero entries from SMA's
sma_v1 = sma_v1(200:length(sma_v1));
sma_v2 = sma_v2(200:length(sma_v2));
sma_v3 = sma_v3(200:length(sma_v3));
sma_v4 = sma_v4(200:length(sma_v4));

% eliminate a bunch of entries from adjusted SMA's
sma_v1_adj = sma_v1_adj(199:length(sma_v1_adj));
sma_v2_adj = sma_v2_adj(199:length(sma_v2_adj));
sma_v3_adj = sma_v3_adj(199:length(sma_v3_adj));
sma_v4_adj = sma_v4_adj(199:length(sma_v4_adj));

% elimate points prior to 10 seconds for flights
noisy_flight1 = noisy_flight1(200:length(noisy_flight1));
noisy_flight2 = noisy_flight2(200:length(noisy_flight2));
noisy_flight3 = noisy_flight3(200:length(noisy_flight3));
noisy_flight4 = noisy_flight4(200:length(noisy_flight4));


R = 287; % ideal gas constant of air
g0 = 9.8; % gravitational constant
a = -0.0065; 
rho_ref = 1.225; 
T_ref = 288.15; 


% TIME FOR SOME DYNAMICS
t = linspace(10,28,1800);
%t = linspace(18.2,28,2000);
URRG = 728; % URRG start in ft
spaceport = 4595; % spaceport start in ft

starting_data_point = 1;

% dry masses from launches
mtl22 = 91.3;
mcl22 = 98.7;
mtl23 = 127.5;
mcl23 = 134.5;

% temperatures from launches
    % TEMP TL 22 = 292.55 K
    % TEMP CL 22 = 304.85 K
    % TEMP TL 23 = 301.25 K
    % TEMP CL 23 = 313.35 K
% adjusted Reference temperatures
T_TL22 = 300 + (292.55-295)*4;
T_CL22 = 300 + (304.85-295)*4;
T_TL23 = 300 + (301.25-295)*4;
T_CL23 = 300 + (313.35-295)*4;


[t,x1] = ode45(@(t,X) pic(t,X,URRG,mtl22,T_TL22),t,[noisy_flight1(starting_data_point+1);sma_v1_adj(starting_data_point)]);
[t2,x2] = ode45(@(t,X) pic(t,X,spaceport,mcl22,T_CL22),t,[noisy_flight2(starting_data_point+1);sma_v2_adj(starting_data_point)]);
[t3,x3] = ode45(@(t,X) pic(t,X,URRG,mtl23,T_TL23),t,[noisy_flight3(starting_data_point+1);sma_v3_adj(starting_data_point)]);
[t4,x4] = ode45(@(t,X) pic(t,X,spaceport,mcl23,T_CL23),t,[noisy_flight4(starting_data_point+1);sma_v4_adj(starting_data_point)]);

% REPETITIVE APOGEE PREDICTION

predictions = [];
predictions2 = [];
for i = 3:length(noisy_flight1)-1
    starting_data_point = i;
    [t,x] = ode45(@(t,X) pic(t,X,URRG,mtl22,T_TL22),t,[noisy_flight1(starting_data_point+1);sma_v1_adj(starting_data_point)]);
    predictions(i) = max(x(:,1));
end

for i = 3:length(noisy_flight2)-1
    starting_data_point = i;
    [t2,x2] = ode45(@(t,X) pic(t,X,spaceport,mcl22,T_CL22),t,[noisy_flight2(starting_data_point+1);sma_v2_adj(starting_data_point)]);
    predictions2(i) = max(x2(:,1));
end

for i = 3:length(noisy_flight3)-1
    starting_data_point = i;
    [t3,x3] = ode45(@(t,X) pic(t,X,URRG,mtl23,T_TL23),t,[noisy_flight3(starting_data_point+1);sma_v3_adj(starting_data_point)]);
    predictions3(i) = max(x3(:,1));
end

for i = 3:length(noisy_flight4)-1
    starting_data_point = i;
    [t4,x4] = ode45(@(t,X) pic(t,X,spaceport,mcl23,T_CL23),t,[noisy_flight4(starting_data_point+1);sma_v4_adj(starting_data_point)]);
    predictions4(i) = max(x4(:,1));
end



predictions = predictions';
predictions2 = predictions2';
predictions3 = predictions3';
predictions4 = predictions4';


predictions_residuals = predictions - 10707;
predictions2_residuals = predictions2 - 10350;
predictions3_residuals = predictions3 - 12315;
predictions4_residuals = predictions4 - 10067;

count = length([predictions;predictions2;predictions3;predictions4]);
good_count = 0;

for i = 1:length(predictions)
    if predictions_residuals(i) < 300
        good_count = good_count + 1;
    end
end

for i = 1:length(predictions2)
    if predictions2_residuals(i) < 300
        good_count = good_count + 1;
    end
end

for i = 1:length(predictions3)
    if predictions3_residuals(i) < 300
        good_count = good_count + 1;
    end
end

for i = 1:length(predictions4)
    if predictions4_residuals(i) < 300
        good_count = good_count + 1;
    end
end

good_percent = good_count / count;
disp(good_percent)

figure
scatter(t(3:1799),predictions_residuals(3:1799),'b');
hold on
scatter(t2(3:1799),predictions2_residuals(3:1799),'r');
scatter(t3(3:1799),predictions3_residuals(3:1799),'g');
scatter(t4(3:1799),predictions4_residuals(3:1799),'c');
xlabel('Time into Flight (s)','FontSize',14)
ylabel('Predicted Apogee - Actual Apogee (ft)','FontSize',14)
title('Residuals of Predictions','FontSize',20)
legend('TL 22 Residuals','CL 22 Residuals','TL 23 Residuals','CL 23 Residuals','FontSize',14)
hold off

function Xdot = pic(t,X,xstart,m,T_ref)
    x = X(1);
    xdot = X(2);
    R = 287;
    %R = 85;
    g0 = 9.81;
    g = 32.17; % ft/s
    a = -0.0065;
    rho_ref = 1.225;
    %rho_ref = 1.6;
    T_ref = 288.15;
    %T_ref = 400;
    Cd = 0.536; % coefficient of drag as predicted on open rocket
    %Cd = 0.77;
    S = pi*0.25^2; % cross sectional area of rocket in ft^2
    %S = pi*0.25^2 + (0.05*S); % cross sectional area of rocket in ft^2

    % MASSES
    % DRY MASS TL 22 = 91.3 lbm
    % DRY MASS CL 22 = 98.7 lbm
    % DRY MASS TL 23 = 127.5 lbm
    % DRY MASS CL 23 = 134.5 lbm

    % TEMPERATURES
    % TEMP TL 22 = 292.55 K
    % TEMP CL 22 = 304.85 K
    % TEMP TL 23 = 301.25 K
    % TEMP CL 23 = 313.35 K



    beta = 0.5*Cd*S;
   

    rho = (1/16.0185)*(rho_ref*(1+(a*((x+xstart)*0.3048)/T_ref)).^((-g0/(a*R))-1)); % ft/lb^3

    altitude = x + xstart;

    %if xstart < 2000
    %    rho = (x+xstart)*(1.1350e-5) - 0.0274;
    %else
    %    rho = 0.033;
    %end
    
    
    %T = 80 - 0.00649*(x+xstart);
    %p = 101.29*((T+273.1)/288.08)^5.256;
    %rho = p / (0.2869*(T+273.1));
    %rho = rho*(1/16.0185);


    xdoubledot = -g - beta/m*rho*xdot^2;
    Xdot = [xdot; xdoubledot];
    
end


