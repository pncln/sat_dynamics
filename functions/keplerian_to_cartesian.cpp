#include "keplerian_to_cartesian.h"
#include "../Constants.h"
#include <vector>

std::vector<double> keplerian_to_cartesian(double a, double e, double i, double Omega, double omega, double nu) {
    double p = a * (1 - e * e);
    double r = p / (1 + e * cos(nu));

    double x = r * (cos(Omega) * cos(omega + nu) - sin(Omega) * sin(omega + nu) * cos(i));
    double y = r * (sin(Omega) * cos(omega + nu) + cos(Omega) * sin(omega + nu) * cos(i));
    double z = r * sin(omega + nu) * sin(i);

    double h = sqrt(G * M_earth * p);
    double vx = (x * h * e / (r * p)) * sin(nu) - (h / r) * (cos(Omega) * sin(omega + nu) + sin(Omega) * cos(omega + nu) * cos(i));
    double vy = (y * h * e / (r * p)) * sin(nu) - (h / r) * (sin(Omega) * sin(omega + nu) - cos(Omega) * cos(omega + nu) * cos(i));
    double vz = (z * h * e / (r * p)) * sin(nu) + (h / r) * cos(omega + nu) * sin(i);

    return {x, y, z, vx, vy, vz};
}
