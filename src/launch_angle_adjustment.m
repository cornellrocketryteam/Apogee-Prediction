clc
clear
close all

x_pos = [0];
y_pos = [6787.4];

vy1 = 654.83;
vx1 = tan(deg2rad(10))*vy1;

x_vel = [vx1];
y_vel = [vy1];
t = [0];

dt = 0.01;

for i = 1:(30*(1/dt))

    [x_next,y_next,vx_next,vy_next,] = propogate(x_pos(i),y_pos(i),x_vel(i),y_vel(i),'URRG',dt);

    x_pos(i+1) = x_next;
    y_pos(i+1) = y_next;

    x_vel(i+1) = vx_next;
    y_vel(i+1) = vy_next;

    t(i+1) = i*dt;

end

apogee = max(y_pos);
disp(apogee)

x_pos = x_pos';
y_pos = y_pos';
x_vel = x_vel';
y_vel = y_vel';
t = t';

plot(x_pos,y_pos,'LineWidth',1.5)
xlabel('x Position (ft)','FontSize',16)
ylabel('y Position (ft)','FontSize',16)
title('Launch Vehicle Trajectory Approximation for TL23','FontSize',16)




function [x_next,y_next,vx_next,vy_next] = propogate(x,y,vx,vy,ystart,dt)

if ystart == 'URRG'
    ystart = 728;
else
    ystart = 4595;
end


g0 = 9.81;
g = 32.17; % ft/s
a = -0.0065;
rho_ref = 1.225;
m = 127.5;
Cd = 0.536;
T_ref = 288.15;
R = 287;

rho = (1/16.0185)*(rho_ref*(1+(a*((y+ystart)*0.3048)/T_ref)).^((-g0/(a*R))-1));
A = pi*0.25^2;
dt = 0.01;


v_mag = sqrt(vx^2 + vy^2);
ax = -(1/(2*m))*Cd*rho*(vx^2+vy^2)*A*vx/v_mag;
ay = -g - (1/(2*m))*Cd*rho*(vx^2+vy^2)*A*vy/v_mag;

x_next = x + vx*dt;
y_next = y + vy*dt;

vx_next = vx + ax*dt;
vy_next = vy + ay*dt;


end