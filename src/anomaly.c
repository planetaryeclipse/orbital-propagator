#include "anomaly.h"
#include <math.h>

#include "stumpff.h"

orbit_t calc_semimajor_axis(state_vec_t sv, orbit_t grav_param) {
    orbit_t r_norm = sqrt(sv.rx * sv.rx + sv.ry * sv.ry + sv.rz * sv.rz);
    orbit_t v_norm = sqrt(sv.vx * sv.vx + sv.vy * sv.vy * sv.vz * sv.vz);

    return 1.0 / (2.0 / r_norm - v_norm * v_norm / grav_param);
}

orbit_t calc_universal_anomaly(state_vec_t sv0, orbit_t grav_param, orbit_t dt) {
    // Calculates the semimajor axis and its reciprocal
    orbit_t a = calc_semimajor_axis(sv0, grav_param);
    orbit_t a_recip = 1.0 / a;

    // Calculates an initial estimate of universal anomaly
    orbit_t x = sqrt(grav_param) * (a_recip < 0 ? -a_recip : a_recip) * dt;

    orbit_t r0_norm = sqrt(sv0.rx * sv0.rx + sv0.ry * sv0.ry + sv0.rz * sv0.rz);

    // Calculates the initial radial velocity
    orbit_t v0_radial_norm = (sv0.rx * sv0.vx + sv0.ry * sv0.vy + sv0.rz * sv0.vz) / r0_norm;

    anomaly_converge_n_t i = 1;
    while (i <= ANOMALY_CONVERGENCE_MAX_N) {
        // To be fed to the Stumpff functions
        stumpff_t z = a_recip * x * x;

        // Calculates a function and derivative derived from the universal formulation for use
        // in Newton-Raphson root finding to convergence at the true value of the universal anomaly

        orbit_t f_i = r0_norm * v0_radial_norm / sqrt(grav_param) * x * x * stumpff_c(z) +
                      (1 - a_recip * r0_norm) * x * x * x *
                      stumpff_s(z) + r0_norm * x -
                      sqrt(grav_param) * dt;

        orbit_t df_i = r0_norm * v0_radial_norm / sqrt(grav_param) * x * (1 - a_recip * x * x * stumpff_s(z)) +
                       (1 - a_recip * r0_norm) * x * x *
                       stumpff_c(z) + r0_norm;

        // Checks if the new change is significant enough to indicate that the convergence of the
        // universal anomaly requires at least one more iteration
        orbit_t ratio = f_i / df_i;
        if (((ratio < 0) ? -ratio : ratio) <= ANOMALY_CONVERGENCE_CRITERION)
            break;

        // Applies the change calculated in this iteration
        x -= ratio;

        ++i;
    }

    // Check for non-convergence with i

    return x;
}

state_vec_t propagate_state_vec(state_vec_t sv0, orbit_t grav_param, orbit_t dt) {
    // Calculates the universal anomaly
    orbit_t x = calc_universal_anomaly(sv0, grav_param, dt);

    orbit_t r0_norm = sqrt(sv0.rx * sv0.rx + sv0.ry * sv0.ry + sv0.rz * sv0.rz);
    orbit_t v0_norm = sqrt(sv0.vx * sv0.vx + sv0.vy * sv0.vy * sv0.vz * sv0.vz);

    // Calculates the semimajor axis and its reciprocal
    orbit_t a = calc_semimajor_axis(sv0, grav_param);
    orbit_t a_recip = 1.0 / a;

    // To be fed to the Stumpff functions
    stumpff_t z = a_recip * x * x;

    // Calculates the Lagrange coefficients to propagate the position
    state_vec_t svt;

    orbit_t lagrange_f = 1.0 - x * x / r0_norm * stumpff_c(z);
    orbit_t lagrange_g = dt - 1.0 / sqrt(grav_param) * x * x * x * stumpff_s(z);

    // Calculates the new position and caches the magnitude
    svt.rx = lagrange_f * sv0.rx + lagrange_g * sv0.vx;
    svt.ry = lagrange_f * sv0.ry + lagrange_g * sv0.vy;
    svt.rz = lagrange_f * sv0.rz + lagrange_g * sv0.vz;

    orbit_t rt_norm = sqrt(svt.rx * svt.rx + svt.ry * svt.ry + svt.rz * svt.rz);

    // Calculates the Lagrange coefficients to propagate the velocity
    orbit_t lagrange_df = sqrt(grav_param) / (rt_norm * r0_norm) * (a_recip * x * x * x * stumpff_s(z) - x);
    orbit_t lagrange_dg = 1.0 * x * x / rt_norm * stumpff_c(z);

    // Computes the new velocity
    svt.vx = lagrange_df * sv0.rx + lagrange_dg * sv0.vx;
    svt.vy = lagrange_df * sv0.ry + lagrange_dg * sv0.vy;
    svt.vz = lagrange_df * sv0.rz + lagrange_dg * sv0.vz;

    return svt;
}