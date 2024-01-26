altitude_polynomial_regression

time1 = linspace(8,28,(28-8)*100);
time2 = linspace(8,28,(28-8)*100);
% create simulated number of data points equivalent to 100 Hz readings

m_altitude1 = model1(1)*time1.^3 + model1(2)*time1.^2 + model1(3)*time1 + model1(4);
m_altitude2 = model2(1)*time2.^3 + model2(2)*time2.^2 + model2(3)*time2 + model2(4);

noisy_flight1 = zeros(1,length(m_altitude1));
noisy_flight2 = zeros(1,length(m_altitude2));

% simulate noisy data
noise_factor = 0; % higher the number, the noisier the data

for i = 1:length(m_altitude1)
    if rand() > 0.5
        noisy_flight1(i) = m_altitude1(i) + noise_factor*rand();
    else
        noisy_flight1(i) = m_altitude1(i) - noise_factor*rand();
    end
end

for i = 1:length(m_altitude2)
    if rand() > 0.5
        noisy_flight2(i) = m_altitude2(i) + noise_factor*rand();
    else
        noisy_flight2(i) = m_altitude2(i) - noise_factor*rand();
    end
end


% obtain velocity from noisy flight data
rolling_L = 20;
smoothed_altitude = zeros(1,length(noisy_flight1) - rolling_L);
smoothed_altitude2 = zeros(1,length(noisy_flight2) - rolling_L);

for i = 1:length(smoothed_altitude)
    smoothed_altitude(i) = mean(noisy_flight1(i:(i+rolling_L)));
    smoothed_altitude2(i) = mean(noisy_flight2(i:(i+rolling_L)));
end

smoothed_velocity = zeros(1,length(smoothed_altitude)-1);
smoothed_velocity2 = zeros(1,length(smoothed_altitude)-1);

for i = 2:length(smoothed_velocity)
    smoothed_velocity(i-1) = (smoothed_altitude(i) - smoothed_altitude(i-1)) / 0.01;
end

for i = 2:length(smoothed_velocity2)
    smoothed_velocity2(i-1) = (smoothed_altitude2(i) - smoothed_altitude2(i-1)) / 0.01;
end


figure
%plot(time1,noisy_flight1)
plot(time1((rolling_L+1):length(time1)),smoothed_altitude,'LineWidth',3)
hold on
%plot(time1(rolling_L+2:length(time1)),smoothed_velocity)
xlabel('Time into Flight (s)','FontSize',14)
ylabel('Altitude (ft)','FontSize',14)
%title('Simulated Altitude Data (TL 22 Scenario)','FontSize',20)
title('Simulated Altitude Data (TL 22 Scenario) with Model','FontSize',20)


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

% TIME FOR SOME DYNAMICS
t = linspace(8.2,28,2000);
%t = linspace(18.2,28,2000);
URRG = 728; % URRG start in ft
spaceport = 4595; % spaceport start in ft

starting_data_point = 1;

mtl22 = 91.3;
mcl22 = 98.7;



[t,x1] = ode45(@(t,X) pic(t,X,URRG,mtl22),t,[smoothed_altitude(starting_data_point+1);smoothed_velocity(starting_data_point)]);
[t2,x2] = ode45(@(t,X) pic(t,X,spaceport,mcl22),t,[smoothed_altitude2(starting_data_point+1);smoothed_velocity2(starting_data_point)]);
plot(t,x1(:,1),'LineWidth',3)

plot(time2((rolling_L+1):length(time2)),smoothed_altitude2,'LineWidth',3)

plot(t2,x2(:,1),'LineWidth',3)
legend('Simulated Altitude Data TL 22','Dynamics Model Prediction 1 (w/ Noise) at 8.2 seconds','Simulated Altitude Data CL 22','Dynamics Model Prediction 2 (w/ Noise)','Location','Southeast','FontSize',15)

hold off

% REPETITIVE APOGEE PREDICTION

predictions = [];
predictions2 = [];
for i = 1:length(smoothed_altitude)-100
    starting_data_point = i;
    [t,x] = ode45(@(t,X) pic(t,X,URRG,mtl22),t,[smoothed_altitude(starting_data_point+1);smoothed_velocity(starting_data_point)]);
    predictions(i) = max(x(:,1));
end

for i = 1:length(smoothed_altitude2)-100
    starting_data_point = i;
    [t2,x2] = ode45(@(t,X) pic(t,X,URRG,mcl22),t,[smoothed_altitude2(starting_data_point+1);smoothed_velocity2(starting_data_point)]);
    predictions2(i) = max(x2(:,1));
end


predictions = predictions';
predictions2 = predictions2';
figure
predictions_residuals = predictions - 10707;
predictions2_residuals = predictions2 - 10350;

scatter(t(1:1799),predictions_residuals(1:1799),'b');
hold on
scatter(t2(1:1799),predictions2_residuals(1:1799),'r');
xlabel('Time into Flight','FontSize',14)
ylabel('Predicted Apogee - Actual Apogee (ft)','FontSize',14)
title('Residuals of Predictions','FontSize',20)
legend('TL 22 Residuals','CL 22 Residuals')
hold off


function Xdot = pic(t,X,xstart,m)
    x = X(1);
    xdot = X(2);
    R = 287;
    g0 = 9.81;
    g = 32.17; % ft/s
    a = -0.0065;
    rho_ref = 1.225;
    T_ref = 288.15;
    Cd = 0.536; % coefficient of drag as predicted on open rocket
    %Cd = 0.77;
    S = pi*0.25^2; % cross sectional area of rocket in ft^2

    % MASSES
    % DRY MASS TL 22 = 91.3 lbm
    % DRY MASS CL 22 = 98.7 lbm
    % DRY MASS TL 23 = 127.5 lbm
    % DRY MASS CL 23 = 134.5 lbm



    beta = 0.5*Cd*S;
   

    rho = (1/16.0185)*(rho_ref*(1+(a*((x+xstart)*0.3048)/T_ref)).^((-g0/(a*R))-1)); % ft/lb^3

    xdoubledot = -g - beta/m*rho*xdot^2;
    Xdot = [xdot; xdoubledot];
    
end
