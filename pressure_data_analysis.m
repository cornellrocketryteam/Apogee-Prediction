clc
clear
close all

data = readtable('Flight_Data.csv');

time = table2array(data(:,1));
altitude = table2array(data(:,2));
pressure = table2array(data(:,3));
velocity = table2array(data(:,4));
flight = table2cell(data(:,9));


figure
scatter(pressure,altitude,'b')
xlabel('Pressure')
ylabel('Altitude')
title('Altitude v. Pressure')

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
        pressure1(w) = pressure(i);
        w = w + 1;
    elseif char(launch) == 'CL 22'
        time2(x) = time(i);
        altitude2(x) = altitude(i);
        pressure2(x) = pressure(i);
        x = x + 1;
    elseif char(launch) == 'TL 23'
        time3(y) = time(i);
        altitude3(y) = altitude(i);
        pressure3(y) = pressure(i);
        y = y + 1;
    else
        time4(z) = time(i);
        altitude4(z) = altitude(i);
        pressure4(z) = pressure(i);
        z = z + 1;
    end

end

delta_p1 = [];
delta_p2 = [];
delta_p3 = [];
delta_p4 = [];


for i = 1:(length(time1)-1)
    delta_p1(i) = (pressure1(i+1) - pressure1(i))/0.05;
end

for i = 1:(length(time2)-1)
    delta_p2(i) = (pressure2(i+1) - pressure2(i))/0.05;
end

for i = 1:(length(time3)-1)
    delta_p3(i) = (pressure3(i+1) - pressure3(i))/0.05;
end

for i = 1:(length(time4)-1)
    delta_p4(i) = (pressure4(i+1) - pressure4(i))/0.05;
end

figure
scatter(time1(1:length(time1)-1),delta_p1)
xlabel('Time into flight')
ylabel('Change in Pressure')
title('Change in Pressure v. Time into Flight')

