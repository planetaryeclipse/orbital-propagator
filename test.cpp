#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

extern "C" {
#include "stumpff.h"
}

using Catch::Approx;

TEST_CASE("Trigonometric Stumpff C(z) function", "[stumpff]"){
    CHECK(stumpff_c(-15) == Approx(1.5368808205));
    CHECK(stumpff_c(0) == Approx(0.5000000000));
    CHECK(stumpff_c(15) == Approx(0.1162830848));
}

TEST_CASE("Trigonometric Stumpff S(z) function", "[stumpff]"){
    CHECK(stumpff_s(-15) == Approx(0.3470095432));
    CHECK(stumpff_s(0) == Approx(0.1666666667));
    CHECK(stumpff_s(15) == Approx(0.0781634938));
}