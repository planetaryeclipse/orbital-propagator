function [r_new, v_new] = propagate_state_vec(r0, v0, g_param, X, dt, stumpff_n)
    r0_mag = norm(r0);
    v0_mag = norm(v0);

    a = calc_semimajor(r0, v0, g_param);
    a_recip = 1 / a;

    z = a_recip * X^2;

    % Calculates the Lagrange coefficients to propagate the position
    lagrange_f = 1 - X^2 / r0_mag * stumpff_c_trig(z);
    lagrange_g = dt - 1 / sqrt(g_param) * X^3 ...
        * stumpff_s_trig(z);

    % Computes the new position and caches the magnitude
    r_new = lagrange_f * r0 + lagrange_g * v0;
    r_new_mag = norm(r_new);

    % Calculates the Lagrange coefficients to propagate the velocity
    lagrange_df = sqrt(g_param) / (r_new_mag * r0_mag) ...
        * (a_recip * X^3 * stumpff_s_trig(z) - X);
    lagrange_dg = 1 - X^2 / r_new_mag * stumpff_c_trig(z);

    % Computes the new velocity
    v_new = lagrange_df * r0 + lagrange_dg * v0;
end