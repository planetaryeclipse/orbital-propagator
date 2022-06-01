function X = calculate_univ_anomaly(r0, v0, g_param, dt, crit, stumpff_n)
    % Calculates the semimjaor axis and its reciprocol
    a = calc_semimajor(r0, v0, g_param);
    a_recip = 1 / a;
    
    % Calculates an initial estimate of X
    X = sqrt(g_param) * abs(a_recip) * dt;
    
    % Caches this value for performance and clarity
    r0_mag = norm(r0);

    % Calculates the initial radial velocity
    v_r0 = dot(r0, v0) / norm(r0);

    while 1==1
        % To be fed into the Stumpff functions
        z = a_recip * X^2;

        % Calculates the function and derivative used for the
        % Newton-Raphson method solver for X
        f_i = r0_mag * v_r0 / sqrt(g_param) ...
                * X^2 * stumpff_c_trig(z) ...
                + (1 - a_recip * r0_mag) ...
                * X^3 * stumpff_s_trig(z) ...
                + r0_mag * X - sqrt(g_param) * dt;

        df_i = r0_mag * v_r0 / sqrt(g_param) * X ...
            * (1 - a_recip * X^2 * stumpff_s_trig(z)) ...
            + (1 - a_recip * r0_mag) ...
            * X^2 * stumpff_c_trig(z) ...
            + r0_mag;
        
        % Checks to see if the new change is significant enough for the
        % calculation of X to use another iteration
        ratio = f_i / df_i;
        if abs(ratio) <= crit
            break;
        end

        % Iterates the next value of X if the process has not terminated
        % and convergence is to continue
        X = X - ratio;
    end
end