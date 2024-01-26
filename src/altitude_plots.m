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

figure
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

scatter(time1,altitude1,'g')
hold on
scatter(time2,altitude2,'b')
scatter(time3,altitude3,'r')
scatter(time4,altitude4,'m')

xlabel('Time (s)','FontSize',16)
ylabel('Altitude (ft)','FontSize',16)
title('Post-Burnout RRC3 Altitude Data for CRT Flights','FontSize',20)
legend('Test Launch 22','Competition Launch 22','Test Launch 23','Competition Launch 23','Location','southeast','FontSize',16)

hold off

