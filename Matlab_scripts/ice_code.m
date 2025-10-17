% Constants
CR = 9.3;
theta0 = 330;
m = 2.5;
a = 6;
deltatheta = 60;
R = 287;
gamma = 1.35;
Vd = 500e-6;
Vc = (Vd / (CR - 1));
T0 = 300;
P0 = 1e5;

% Initialize arrays
thetadeg_all = 0:2.5:720;
T_all = zeros(size(thetadeg_all));
P_all = zeros(size(thetadeg_all));

% Loop through crank angles
for i = 1:length(thetadeg_all)
thetadeg = thetadeg_all(i);
if thetadeg < 720
theta = thetadeg;
else
theta = mod(thetadeg, 720);
end
if theta < theta0
xb = 0;
elseif theta > (theta0 + deltatheta)
xb = 1;
else

xb = 1 - exp(-a * ((theta - theta0) / deltatheta)^m);
end
Qtotal = 300;
Q = Qtotal * xb;
Vi = Vc + 0.5 * Vd * (1 - cosd(theta0));
Ti = T0 * (((Vc + Vd) / Vi)^(gamma - 1));
Pi = P0 * (((Vc + Vd) / Vi)^gamma);
Vj = Vc + 0.5 * Vd * (1 - cosd(theta0 + deltatheta));
Tj = Ti + (((gamma - 1) * 300) / ((Vc + Vd) * R));
Pj = Pi * (Tj / T0) * ((Vc / Vj)^gamma);
Vk = Vc + 0.5 * Vd * (1 - cosd(540));
if theta <= 180
T = T0;
P = P0;
V = Vc + 0.5 * Vd * (1 - cosd(theta));
elseif (theta > 180 && theta <= theta0)
V = Vc + 0.5 * Vd * (1 - cosd(theta));
T = T0 * (((Vc + Vd) / V)^(gamma - 1));
P = P0 * (((Vc + Vd) / V)^gamma);
elseif (theta > theta0 && theta <= theta0 + deltatheta)
V = Vc + 0.5 * Vd * (1 - cosd(theta));
T = Ti + (((gamma - 1) * Q) / ((Vc+Vd) * R));
P = Pi * (T / T0) * ((Vc / V)^gamma);
elseif (theta0 + deltatheta < theta && theta <= 540)
V = Vc + 0.5 * Vd * (1 - cosd(theta));
T = Tj * ((Vj / V)^(gamma - 1));
P = Pj * ((Vj / V)^gamma);
elseif (theta > 540 && theta <= 720)
T = Tj * ((Vj / Vk)^(gamma - 1));
P = Pj * ((Vj / Vk)^gamma);
V = Vc + 0.5 * Vd * (1 - cosd(theta));
end
% Store the values

T_all(i) = T;
P_all(i) = P;
end
% Plot Temperature vs. Crank Angle
figure;
scatter(thetadeg_all, T_all, 15, 'r', 'filled');
xlabel('Crank Angle (degrees)');
ylabel('Temperature (K)');
title('Temperature vs. Crank Angle');
grid on;
% Plot Pressure vs. Crank Angle
figure;
scatter(thetadeg_all, P_all / 1e5, 15, 'k', 'filled'); % Converted to bar
xlabel('Crank Angle (degrees)');
ylabel('Pressure (bar)');
title('Pressure vs. Crank Angle');
grid on;

% Engine and Combustion Constants
CR = 9.3; % Compression ratio
pulseWidth = 0.72; % Duration of laser pulse in degrees
pulseIntensity = 0.1; % Laser fraction of total heat
theta0 = 330; % Start of combustion
a = 6; % Wiebe parameter (shape)
m = 2.5; % Wiebe exponent
deltatheta = 60; % Duration of combustion
R = 287; % Specific gas constant for air (J/kg·K)
gamma = 1.32; % Specific heat ratio (typical for combustion gases)
Vd = 500e-6; % Displacement volume (m^3)
Vc = Vd / (CR - 1); % Clearance volume (m^3)
T0 = 300; % Initial temperature (K)
P0 = 1e5; % Initial pressure (Pa)
Qtotal = 300; % Total heat released (J)

efficiency = 0.5; % Effective conversion to thermal energy

% Exhaust conditions (post 540 degrees)
T_exhaust = 500; % Exhaust temperature (K)
P_exhaust = 1e5; % Exhaust pressure (Pa) ~ 1 bar

% Arrays for results
thetadeg_all = 0:2.5:720;
T_all = zeros(size(thetadeg_all));
P_all = zeros(size(thetadeg_all));

% Variables to track the pressure and temperature in exhaust phase
T_exhaust_value = NaN; % Initialize a variable for the exhaust temperature
P_exhaust_value = NaN; % Initialize a variable for the exhaust pressure

% Main loop over crank angles
for i = 1:length(thetadeg_all)
thetadeg = thetadeg_all(i);
theta = mod(thetadeg, 720); % Ensure theta is within 0-720 degrees

% Instantaneous volume calculation
V = Vc + 0.5 * Vd * (1 - cosd(theta));

% Initial conditions
Vi = Vc + 0.5 * Vd * (1 - cosd(theta0));
Ti = T0 * (((Vc + Vd) / Vi)^(gamma - 1));
Pi = P0 * (((Vc + Vd) / Vi)^gamma);
if theta <= 180
T = T0;
P = P0;
elseif theta > 180 && theta < theta0
T = T0 * (((Vc + Vd) / V)^(gamma - 1));
P = P0 * (((Vc + Vd) / V)^gamma);

elseif theta >= theta0 && theta <= theta0 + pulseWidth
% Laser-assisted rapid combustion
Q = pulseIntensity * Qtotal * ((theta - theta0) / pulseWidth);
T = Ti + (efficiency * (gamma - 1) * Q) / (V * R);
P = Pi * (T / T0) * ((Vc / V)^gamma);
else
% Gradual combustion via Wiebe function
xb = 1 - exp(-a * ((theta - (theta0 + pulseWidth)) / deltatheta)^m);
xb = min(max(xb, 0), 1);
Q = pulseIntensity * Qtotal + (1 - pulseIntensity) * Qtotal * xb;
T = Ti + (efficiency * (gamma - 1) * Q) / (V * R);
P = Pi * (T / T0) * ((Vc / V)^gamma);
end
% After 540 degrees (exhaust phase)
if theta > 540
% If T_exhaust_value and P_exhaust_value are not initialized yet, store them
if isnan(T_exhaust_value)
T_exhaust_value = T;
end
if isnan(P_exhaust_value)
P_exhaust_value = P;
end
% Keep temperature and pressure constant in exhaust phase
T = T_exhaust_value;
P = P_exhaust_value;
end
% Store results
T_all(i) = T;
P_all(i) = P;
end
% Plot Temperature
figure;
scatter(thetadeg_all, T_all, 15, 'r', 'filled');

xlabel('Crank Angle (degrees)');
ylabel('Temperature (K)');
title('Temperature vs. Crank Angle');
grid on;
% Plot Pressure
figure;
scatter(thetadeg_all, P_all / 1e5, 15, 'k', 'filled'); % Convert from Pa to bar
xlabel('Crank Angle (degrees)');
ylabel('Pressure (bar)');
title('Pressure vs. Crank Angle');
grid on;