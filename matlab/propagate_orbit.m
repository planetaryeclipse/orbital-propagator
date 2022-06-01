% This script tests and displays the result of the propagation of an orbit
% using the universal formulation propagator vs the solving of the
% direct differential equation of the orbit.

% Note that the use of the Stumpff trigonometric equivalent functions are
% used as the infinite series requires far more terms to converge at higher
% values of the universal anomaly. Additionally, at higher timesteps, the
% universal anomaly appears to take longer to converge.

% Physical constants
grav_constant = 6.6743e-11; % m^3 kg^-1 s^-2
earth_radius = 6.378e6; % m
earth_mass = 5.974e24; % kg

% Constants of the 2-body system
g_param = grav_constant * earth_mass; % Satellite mass is negligible

% Initial state vector of the satellite
r0 = [200e3 + earth_radius ; 0; 0];
v0 = [0; 10000; 2000]; % sqrt(GM/r)

% Parameters for calculating the orbit
crit = 1e-8;
stumpff_n = 20; % Need more terms at higher X for proper orbit simulation
timestep = 60; % Timestep between computed positions (every minute)
num_points = 10000; % Number of orbital positions computed

% Calculates positions using the universal anomaly
univ_positions = zeros(3, num_points);
for i=0:num_points-1
    X = calculate_univ_anomaly(r0, v0, g_param, ...
        timestep * i, crit, stumpff_n);
    [r, v] = propagate_state_vec(r0, v0, g_param, X, ...
        timestep * i, stumpff_n);
    
    univ_positions(1, i + 1) = r(1);
    univ_positions(2, i + 1) = r(2);
    univ_positions(3, i + 1) = r(3);

    fprintf("Computing position: %g\n", i);
    % fprintf("\tVector: x: %g, y: %g, z: %g\n", r(1), r(2), r(3));
end

% Calculates positions using the ODE solver

% Represent ODE as a series of linear first order ODEs
y0 = [r0(1) ; v0(1) ; r0(2) ; v0(2) ; r0(3) ; v0(3)];
[t, y] = ode89(@(t,y) orbital_ode_f(t, y, g_param), ...
    [0, num_points * timestep], y0);

% Displays the computed orbits
camproj("perspective");
set(gcf, 'position', [100, 100, 800, 800]);
hold on;
plot3(y(:,1), y(:, 3), y(:, 5), 'g', 'MarkerSize', 10);
plot3(univ_positions(1, :), univ_positions(2, :), ...
    univ_positions(3, :), 'r', 'MarkerSize', 10);
plot3(0:0, 0:0, 0:0, 'ob', 'MarkerSize', 10);

daspect([1 1 1]);
pbaspect([2 1 1]); % Modify this to change the box containing the orbit