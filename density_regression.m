% air density

clear
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
% adjusted reference temperatures

% modify lengths of velocity vectors
sma_v1_adj = sma_v1_adj(2:1682);
sma_v2_adj = sma_v2_adj(2:1606);
sma_v3_adj = sma_v3_adj(2:1800);
sma_v4_adj = sma_v4_adj(2:1645);

alt1 = noisy_flight1(2:length(sma_v1_adj)+1);
alt2 = noisy_flight2(2:length(sma_v2_adj)+1);
alt3 = noisy_flight3(2:length(sma_v3_adj)+1);
alt4 = noisy_flight4(2:length(sma_v4_adj)+1);

Cd = 0.536;
S = pi*0.25^2;

rho1 = [];
rho2 = [];
rho3 = [];
rho4 = [];

rho1 = density(alt1,sma_v1_adj,mtl22,URRG)';
rho2 = density(alt2,sma_v2_adj,mcl22,spaceport)';
rho3 = density(alt3,sma_v3_adj,mtl23,URRG)';
rho4 = density(alt4,sma_v4_adj,mcl23,spaceport)';


x = 800;
figure
plot(alt1(1:length(alt1)-x)+URRG,rho1)
hold on
plot(alt2(1:length(alt2)-x)+spaceport,rho2)
plot(alt3(1:length(alt3)-x)+URRG,rho3)
plot(alt4(1:length(alt4)-x)+spaceport,rho4)

xlabel('Altitude (ft)')
ylabel('Density')

legend('TL 22','CL 22','TL 23','CL 23')

hold off

alt1 = alt1(1:length(alt1)-x)'+URRG;
alt2 = alt2(1:length(alt2)-x)'+spaceport;
alt3 = alt3(1:length(alt3)-x)'+URRG;
alt4 = alt4(1:length(alt4)-x)'+spaceport;

alt_URRG = [alt1;alt3];
rho_URRG = [rho1;rho3];

alt_spaceport = [alt2;alt4];
rho_spaceport = [rho2;rho4];






function rho_vector = density(alt_vector,v_vector,mass,location)
    Cd = 0.536;
    S = pi*0.25^2;
    g = 32.17;
    t = 0.01;
    m = mass;
    rho_vector = [];

    disp('NEW FLIGHT')

    for i = 1:length(alt_vector)-800
        y0 = alt_vector(i) + location;
        y = alt_vector(i+1) + location;
        v0 = v_vector(i);
        v = v_vector(i+1);

        a = (v^2 - v0^2) / (2*(y-y0));
        disp(a)

        rho_vector(i) = -2*m*(a+g) / (Cd*S*v^2);
        



    end

end
    





