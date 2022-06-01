#ifndef UNIVERSAL_PROPAGATOR_STUMPFF_H
#define UNIVERSAL_PROPAGATOR_STUMPFF_H

#define STUMPFF_MAX_VALID_N 12

typedef double stumpff_t;
typedef unsigned char stumpff_n_t;
typedef unsigned long stumpff_factorial_t;

stumpff_factorial_t stumpff_factorial(stumpff_factorial_t n);

stumpff_t stumpff_c(stumpff_t z);

stumpff_t stumpff_s(stumpff_t z);

stumpff_t stumpff_c_series(stumpff_t z, stumpff_n_t n);

stumpff_t stumpff_s_series(stumpff_t z, stumpff_n_t n);

#endif //UNIVERSAL_PROPAGATOR_STUMPFF_H
