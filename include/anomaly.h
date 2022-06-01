#ifndef UNIVERSAL_PROPAGATOR_ANOMALY_H
#define UNIVERSAL_PROPAGATOR_ANOMALY_H

#define ANOMALY_CONVERGENCE_CRITERION 1e-8
#define ANOMALY_CONVERGENCE_MAX_N 10

typedef double orbit_t;
typedef double anomaly_t;

typedef unsigned char anomaly_converge_n_t;

typedef struct {
    orbit_t rx, ry, rz;
    orbit_t vx, vy, vz;
} state_vec_t;

orbit_t calc_semimajor_axis(state_vec_t sv, orbit_t grav_param);

orbit_t calc_universal_anomaly(state_vec_t sv0, orbit_t grav_param, orbit_t dt);

state_vec_t propagate_state_vec(state_vec_t sv0, orbit_t grav_param, orbit_t dt);

#endif //UNIVERSAL_PROPAGATOR_ANOMALY_H
