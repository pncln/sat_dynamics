#include "cartesian_to_keplerian.h"
#include <cmath>
#include "../Constants.h"
#include "cross_product.h"
#include "magnitude.h"
#include "dot_product.h"
#include "scale_vector.h"
#include "subtract_vectors.h"

std::vector<double> cartesian_to_keplerian(double x, double y, double z, double vx, double vy, double vz) {
    std::vector<double> r = {x, y, z};
    std::vector<double> v = {vx, vy, vz};

    double r_mag = sqrt(x*x + y*y + z*z);
    double v_mag = sqrt(vx*vx + vy*vy + vz*vz);

    std::vector<double> h = cross_product(r, v);
    double h_mag = magnitude(h);

    double e_scalar = (v_mag*v_mag - G*M_earth/r_mag) * r_mag - dot_product(r, v) * dot_product(r, v);
    e_scalar /= G * M_earth;
    std::vector<double> e_vec = subtract_vectors(scale_vector(subtract_vectors(scale_vector(v, r_mag), scale_vector(r, dot_product(r, v)/r_mag)), 1.0/G/M_earth), scale_vector(r, 1.0/r_mag));
    double e = magnitude(e_vec);

    double a = -G * M_earth / (2 * (v_mag*v_mag/2 - G*M_earth/r_mag));
    double i = acos(h[2] / h_mag);
    double Omega = atan2(h[0], -h[1]);
    double omega = atan2(e_vec[2]/sin(i), e_vec[0]*cos(Omega) + e_vec[1]*sin(Omega));
    double nu = atan2(dot_product(r, v) * h_mag, h_mag*h_mag - r_mag*G*M_earth);

    return {a, e, i, Omega, omega, nu};
}
