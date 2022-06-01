function dydt = orbital_ode_f(t, y, g_param)
    % Distance of satellite needed for nonlinear ODE

    r = sqrt(y(1)^2 + y(3)^2 + y(5)^2);

    dydt = zeros(6, 1);
    dydt(1) = y(2);
    dydt(2) = - g_param / r^3 * y(1);
    dydt(3) = y(4);
    dydt(4) = - g_param / r^3 * y(3);
    dydt(5) = y(6);
    dydt(6) = - g_param / r^3 * y(5);
end