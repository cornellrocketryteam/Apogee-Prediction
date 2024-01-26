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

TL_model = polyfit(TL_rrc3v,TL_d,3);
CL_model = polyfit(CL_rrc3v,CL_d,3);

TLMP = TL_model(1)*TL_rrc3v.^3 + TL_model(2)*TL_rrc3v.^2 + TL_model(3)*TL_rrc3v + TL_model(4); % test launch model predictions based on RRC3 velocities
CLMP = CL_model(1)*CL_rrc3v.^3 + CL_model(2)*CL_rrc3v.^2 + CL_model(3)*CL_rrc3v + CL_model(4);

figure % figure 1

scatter(rrc3_v1,d1,'g')
hold on

scatter(rrc3_v2,d2,'r')

scatter(rrc3_v3,d3,'b')

scatter(rrc3_v4,d4,'c')

hold off

TL_residuals = TL_d - TLMP;
CL_residuals = CL_d - CLMP;

xlabel('RRC3 Filtered Velocity','FontSize',16)
ylabel('Distance to Apogee','FontSize',16)
title('Post-Burnout RRC3 Altitude Data for CRT Flights','FontSize',20)
legend('Test Launch 22','Competition Launch 22','Test Launch 23','Competition Launch 23','Location','southeast','FontSize',16)

figure
scatter(TL_rrc3v,TL_residuals)
xlabel('RRC3 Filtered Velocity (ft/s)')
ylabel('Actual Distance to Apogee - Predicted Distance to Apogee')
title('Test Launch Apogee Prediction Residuals v. Filtered Velocity','FontSize',18)

figure
scatter(CL_rrc3v,CL_residuals)
xlabel('RRC3 Filtered Velocity (ft/s)')
ylabel('Actual Distance to Apogee - Predicted Distance to Apogee')
title('Competition Launch Apogee Prediction Residuals v. Filtered Velocity','FontSize',18)





