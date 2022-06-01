% Orbital solver with ODE

% r''(t) = - mu/r^3 r(t)

% Physical constants
grav_constant = 6.6743e-11; % m^3 kg^-1 s^-2
earth_radius = 6.378e6; % m
earth_mass = 5.974e24; % kg

% Constants of the 2-body system
g_param = grav_constant * earth_mass; % Satellite mass is negligible

% Initial state vector of the satellite
r0 = [200e3 + earth_radius ; 0; 0];
v0 = [0; 10000*1.05; 2000]; % sqrt(GM/r)

% Represent ODE as a series of linear first order ODEs
y0 = [r0(1) ; v0(1) ; r0(2) ; v0(2) ; r0(3) ; v0(3)];

[t, y] = ode89(@(t,y) orbital_ode_f(t, y, g_param), [0, 100000*60], y0);

% Displays the computed orbits
camproj("perspective");
set(gcf, 'position', [100, 100, 800, 800]);
plot3(y(:,1), y(:, 3), y(:, 5), 'r', 'MarkerSize', 10);
hold on;
plot3(0:0, 0:0, 0:0, 'ob', 'MarkerSize', 10);

daspect([1 1 1]);
pbaspect([1 1 1]);