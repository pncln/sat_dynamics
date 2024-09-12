#ifndef SATELLITE_MODEL_H
#define SATELLITE_MODEL_H

#include <vector>
#include "Constants.h"

class SatelliteModel {
private:
    std::vector<double> state;
    double mass;
    std::vector<double> magnetic_moment;
    std::vector<double> desired_state;

public:
    SatelliteModel(const std::vector<double>& initial_state, double satellite_mass, const std::vector<double>& initial_magnetic_moment);

    void update(double dt, const std::vector<double>& torque, const std::vector<double>& external_disturbance, double t);
    std::vector<double> calculate_derivatives(const std::vector<double>& state, const std::vector<double>& torque, const std::vector<double>& external_disturbance, double t);
    // std::vector<double> keplerian_to_cartesian(double a, double e, double i, double Omega, double omega, double nu);
    std::vector<double> cartesian_to_keplerian(double x, double y, double z, double vx, double vy, double vz);
    std::vector<double> calculate_lunar_perturbation(double a, double e, double i, double Omega, double omega, double nu, double t);
    std::vector<double> calculate_solar_planetary_perturbations(double a, double e, double i, double Omega, double omega, double nu, double t);
    std::vector<double> cartesian_to_keplerian_rates(double a, double e, double i, double Omega, double omega, double nu, const std::vector<double>& accel);
    std::vector<double> calculate_gravitational_torque(double a, double e, double i, double Omega, double omega, double nu);
    std::vector<double> calculate_magnetic_field(double a, double e, double i, double Omega, double omega, double nu, double t);
    std::vector<double> calculate_gravitational_harmonics(double a, double e, double i, double Omega, double omega, double nu);
    std::vector<double> calculate_relativistic_corrections(double a, double e, double i, double Omega, double omega, double nu);
    std::vector<double> calculate_gravitational_force(double x, double y, double z);
    void calculate_orbit_state_vectors(double a, double e, double i, double Omega, double omega, double nu, std::vector<double>& position, std::vector<double>& velocity);
    std::vector<std::vector<double>> matrix_multiply(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B);
    double legendre_polynomial(int n, int m, double x);
    double legendre_polynomial_derivative(int n, int m, double x);

    std::vector<double> generate_elliptical_orbit();

    std::vector<double> matrix_vector_multiply(const std::vector<std::vector<double>>& A, const std::vector<double>& v);

    const std::vector<double>& get_state() const;
    void set_desired_state(const std::vector<double>& new_desired_state);
};

#endif // SATELLITE_MODEL_H
