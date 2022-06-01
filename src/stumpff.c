#include "stumpff.h"
#include <math.h>
#include <stdbool.h>

inline stumpff_factorial_t stumpff_factorial(stumpff_factorial_t n) {
    stumpff_factorial_t val = 1;
    for (int i = 2; i <= n; i++)
        val *= i;

    return val;
}

inline stumpff_t stumpff_c(stumpff_t z) {
    if (z < 0)
        return (cosh(sqrt(-z)) - 1) / (-z);
    else if (z == 0)
        return 1.0 / 2;

    // z > 0
    return (1 - cos(sqrt(z))) / z;
}

inline stumpff_t stumpff_s(stumpff_t z) {
    if (z < 0) {
        stumpff_t sqrt_neg_z = sqrt(-z);
        return (sinh(sqrt_neg_z) - sqrt_neg_z) / (sqrt_neg_z * sqrt_neg_z * sqrt_neg_z);
    } else if (z == 0)
        return 1.0 / 6;

    // z > 0
    stumpff_t sqrt_z = sqrt(z);
    return (sqrt_z - sin(sqrt_z)) / (sqrt_z * sqrt_z * sqrt_z);
}

inline stumpff_t stumpff_c_series(stumpff_t z, stumpff_n_t n) {
    stumpff_t sum = 0;

    stumpff_n_t valid_n = (n < STUMPFF_MAX_VALID_N) ? n : STUMPFF_MAX_VALID_N;
    for (stumpff_n_t i = 0; i < valid_n; i++) {
        bool neg_factor = (i % 2 == 0);
        sum += (neg_factor ? -1 : 1) * pow(z, i) / (stumpff_t) stumpff_factorial(2 * i + 2);
    }

    return sum;
}

inline stumpff_t stumpff_s_series(stumpff_t z, stumpff_n_t n) {
    stumpff_t sum = 0;

    stumpff_n_t valid_n = (n < STUMPFF_MAX_VALID_N) ? n : STUMPFF_MAX_VALID_N;
    for (stumpff_n_t i = 0; i < valid_n; i++) {
        bool neg_factor = (i % 2 == 0);
        sum += (neg_factor ? -1 : 1) * pow(z, i) / (stumpff_t) stumpff_factorial(2 * i + 3);
    }

    return sum;
}