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
        w = w + 1;
    elseif char(launch) == 'CL 22'
        time2(x) = time(i);
        altitude2(x) = altitude(i);
        x = x + 1;
    elseif char(launch) == 'TL 23'
        time3(y) = time(i);
        altitude3(y) = altitude(i);
        y = y + 1;
    else
        time4(z) = time(i);
        altitude4(z) = altitude(i);
        z = z + 1;
    end

end

v1 = [];
v2 = [];
v3 = [];
v4 = [];


for i = 1:(length(time1)-1)
    v1(i) = (altitude1(i+1) - altitude1(i))/0.05;
end

for i = 1:(length(time2)-1)
    v2(i) = (altitude2(i+1) - altitude2(i))/0.05;
end

for i = 1:(length(time3)-1)
    v3(i) = (altitude3(i+1) - altitude3(i))/0.05;
end

for i = 1:(length(time4)-1)
    v4(i) = (altitude4(i+1) - altitude4(i))/0.05;
end

figure
scatter(time1(1:length(time1)-1),v1,'b')
%hold on
%scatter(time2(1:length(time2)-1),v2,'g')
%scatter(time3(1:length(time3)-1),v3,'b')
%scatter(time4(1:length(time4)-1),v4,'m')

xlabel('Time into Flight (s)','FontSize',16)
ylabel('Velocity (ft/s)','FontSize',16)
title('Velocity v. Time into Flight for Test Launch 22','FontSize',20)
%legend('Test Launch 22','Competition Launch 22','Test Launch 23','Competition Launch 23','Location','northeast','FontSize',12)