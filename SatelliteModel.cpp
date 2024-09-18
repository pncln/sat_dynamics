#include <iostream>
#include <cassert>
#include "SatelliteModel.h"
#include "functions/keplerian_to_cartesian.h"
#include "functions/cartesian_to_keplerian.h"
#include "functions/cross_product.h"
#include "functions/magnitude.h"
#include "functions/dot_product.h"
#include "functions/subtract_vectors.h"
#include "functions/add_vectors.h"
#include "functions/multiply_vector.h"
#include "functions/scale_vector.h"


    SatelliteModel::SatelliteModel(const std::vector<double>& initial_state, double satellite_mass, const std::vector<double>& initial_magnetic_moment)
    : state(initial_state), mass(satellite_mass), magnetic_moment(initial_magnetic_moment) {
    // Ensure the initial state vector has the correct size
        assert(state.size() == 13);
    }

    void SatelliteModel::set_desired_state(const std::vector<double>& new_desired_state) {
        desired_state = new_desired_state;
    }

    void SatelliteModel::update(double dt, const std::vector<double>& torque, const std::vector<double>& external_disturbance, double t) {
        std::vector<double> k1 = calculate_derivatives(state, torque, external_disturbance, t);
        std::vector<double> k2 = calculate_derivatives(add_vectors(state, multiply_vector(k1, dt/2)), torque, external_disturbance, t + dt/2);
        std::vector<double> k3 = calculate_derivatives(add_vectors(state, multiply_vector(k2, dt/2)), torque, external_disturbance, t + dt/2);
        std::vector<double> k4 = calculate_derivatives(add_vectors(state, multiply_vector(k3, dt)), torque, external_disturbance, t + dt);

        state = add_vectors(state, multiply_vector(add_vectors(k1, multiply_vector(add_vectors(k2, k3), 2), k4), dt/6));

        // Normalize quaternion
        double norm = sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5] + state[6]*state[6]);
        state[3] /= norm;
        state[4] /= norm;
        state[5] /= norm;
        state[6] /= norm;

        // Ensure Keplerian elements remain within valid ranges
        state[8] = std::max(0.0, std::min(state[8], 0.99999)); // Eccentricity
        state[9] = std::fmod(state[9], 2 * M_PI); // Inclination
        state[10] = std::fmod(state[10], 2 * M_PI); // RAAN
        state[11] = std::fmod(state[11], 2 * M_PI); // Argument of periapsis
        state[12] = std::fmod(state[12], 2 * M_PI); // True anomaly
    }


    const std::vector<double>& SatelliteModel::get_state() const {
        return state;
    }

    std::vector<double> SatelliteModel::calculate_derivatives(const std::vector<double>& state, const std::vector<double>& torque, const std::vector<double>& external_disturbance, double t) {
        std::vector<double> derivatives(13);

        double wx = state[0], wy = state[1], wz = state[2];
        double qx = state[3], qy = state[4], qz = state[5], q0 = state[6];
        double a = state[7], e = state[8], i = state[9], Omega = state[10], omega = state[11], nu = state[12];

        // Angular acceleration (unchanged)
        derivatives[0] = (torque[0] + external_disturbance[0] + (Iy - Iz) * wy * wz) / Ix;
        derivatives[1] = (torque[1] + external_disturbance[1] + (Iz - Ix) * wz * wx) / Iy;
        derivatives[2] = (torque[2] + external_disturbance[2] + (Ix - Iy) * wx * wy) / Iz;

        //  std::cout << "Derivatives calculated: " << derivatives[0] << ", " << derivatives[1] << ", " << derivatives[2] << "..." << std::endl;

        // Quaternion derivatives (unchanged)
        derivatives[3] = 0.5 * (q0 * wx - qz * wy + qy * wz);
        derivatives[4] = 0.5 * (qz * wx + q0 * wy - qx * wz);
        derivatives[5] = 0.5 * (-qy * wx + qx * wy + q0 * wz);
        derivatives[6] = -0.5 * (qx * wx + qy * wy + qz * wz);

        // Keplerian element derivatives
        double n = std::sqrt(G * M_earth / (a * a * a)); // Mean motion
        double p = a * (1 - e * e);
        double h = std::sqrt(G * M_earth * p);

        derivatives[7] = 0; // Semi-major axis rate (assuming no thrust)
        derivatives[8] = 0; // Eccentricity rate (assuming no perturbations)
        derivatives[9] = 0; // Inclination rate (assuming no out-of-plane forces)
        derivatives[10] = 0; // RAAN rate (assuming no out-of-plane forces)
        derivatives[11] = 0; // Argument of periapsis rate (assuming no perturbations)
        derivatives[12] = n * std::pow((1 + e * std::cos(nu)) / (1 - e * e), 2); // True anomaly rate

        // Add perturbations (simplified for demonstration)
        std::vector<double> cartesian = keplerian_to_cartesian(a, e, i, Omega, omega, nu);
        std::vector<double> lunar_pert = calculate_lunar_perturbation(a, e, i, Omega, omega, nu, t);
        std::vector<double> solar_planetary_pert = calculate_solar_planetary_perturbations(a, e, i, Omega, omega, nu, t);


        // Convert perturbations to Keplerian rates (simplified)
        double perturb_magnitude = magnitude(lunar_pert) + magnitude(solar_planetary_pert);
        derivatives[7] += perturb_magnitude * 1e-6; // Approximate effect on semi-major axis
        derivatives[8] += perturb_magnitude * 1e-8; // Approximate effect on eccentricity
        derivatives[11] += perturb_magnitude * 1e-7; // Approximate effect on argument of periapsis

        return derivatives;
    }

    std::vector<double> calculate_quaternion_derivatives(double qx, double qy, double qz, double q0, double wx, double wy, double wz) {
        double qx_dot = 0.5 * (q0 * wx - qz * wy + qy * wz);
        double qy_dot = 0.5 * (qz * wx + q0 * wy - qx * wz);
        double qz_dot = 0.5 * (-qy * wx + qx * wy + q0 * wz);
        double q0_dot = -0.5 * (qx * wx + qy * wy + qz * wz);
        return {qx_dot, qy_dot, qz_dot, q0_dot};
    }

    std::vector<double> SatelliteModel::calculate_gravitational_torque(double a, double e, double i, double Omega, double omega, double nu) {
        std::vector<double> pos = keplerian_to_cartesian(a, e, i, Omega, omega, nu);
        double x = pos[0], y = pos[1], z = pos[2];
        double r = std::sqrt(x*x + y*y + z*z);
        double F_grav = G * M_earth * mass / (r * r);
        double Tx = 3 * F_grav * (y*z) * (Iz - Iy) / (r*r);
        double Ty = 3 * F_grav * (z*x) * (Ix - Iz) / (r*r);
        double Tz = 3 * F_grav * (x*y) * (Iy - Ix) / (r*r);
        return {Tx, Ty, Tz};
    }

    std::vector<double> calculate_gyroscopic_torque(double wx, double wy, double wz) {
        double Tx = (Iy - Iz) * wy * wz;
        double Ty = (Iz - Ix) * wz * wx;
        double Tz = (Ix - Iy) * wx * wy;
        return {Tx, Ty, Tz};
    }

    std::vector<double> SatelliteModel::calculate_gravitational_force(double x, double y, double z) {
        double r = std::sqrt(x*x + y*y + z*z);
        const double epsilon = 1e-10; // Small value to prevent division by zero
        if (r < epsilon) {
            std::cout << "Error: Division by zero in calculate_gravitational_force" << std::endl;
            return {0.0, 0.0, 0.0}; // Return zero force if r is too small
        }
        double F = -G * M_earth * mass / (r * r);
        return {F * x / r, F * y / r, F * z / r};
    }

    std::vector<double> SatelliteModel::calculate_gravitational_harmonics(double a, double e, double i, double Omega, double omega, double nu) {
        std::vector<double> pos = keplerian_to_cartesian(a, e, i, Omega, omega, nu);
        double x = pos[0], y = pos[1], z = pos[2];

        const double J2 = -1.08262617352e-3;  // J2 coefficient
        const double J3 = -2.53265648533e-6;  // J3 coefficient
        const double J4 = -1.61962159137e-6;  // J4 coefficient

        double r = std::sqrt(x*x + y*y + z*z);
        double r_squared = r * r;
        double r_cubed = r * r_squared;
        double r_fourth = r_squared * r_squared;
        double r_fifth = r * r_fourth;
        double r_sixth = r_squared * r_fourth;
        double r_seventh = r * r_sixth;

        double coeff_j2 = -1.5 * J2 * G * M_earth * R_earth * R_earth / r_fifth;
        double coeff_j3 = -2.5 * J3 * G * M_earth * R_earth * R_earth * R_earth / r_seventh;
        double coeff_j4 = -1.875 * J4 * G * M_earth * R_earth * R_earth * R_earth * R_earth / (r_seventh * r);

        double z_squared = z * z;
        double z_cubed = z * z_squared;
        double z_fourth = z_squared * z_squared;

        double x_pert = x * (coeff_j2 * (5 * z_squared / r_squared - 1) +
                             coeff_j3 * (7 * z_cubed / r_cubed - 3 * z / r) * (R_earth / r) +
                             coeff_j4 * (3 - 42 * z_squared / r_squared + 63 * z_fourth / r_fourth));

        double y_pert = y * (coeff_j2 * (5 * z_squared / r_squared - 1) +
                             coeff_j3 * (7 * z_cubed / r_cubed - 3 * z / r) * (R_earth / r) +
                             coeff_j4 * (3 - 42 * z_squared / r_squared + 63 * z_fourth / r_fourth));

        double z_pert = z * (coeff_j2 * (5 * z_squared / r_squared - 3) +
                             coeff_j3 * (7 * z_cubed / r_cubed - 5 * z / r) * (R_earth / r) +
                             coeff_j4 * (15 - 70 * z_squared / r_squared + 63 * z_fourth / r_fourth));

        return {x_pert, y_pert, z_pert};
    }

    std::vector<double> SatelliteModel::generate_elliptical_orbit() {
        double semi_major_axis = R_earth + 10000000; // 10,000 km above Earth's surface
        double eccentricity = 0.1; // Elliptical orbit
        double inclination = 0.5; // 28.6 degrees
        double argument_of_perigee = 0.3; // 17.2 degrees
        double right_ascension = 0.7; // 40.1 degrees
        double true_anomaly = 0.0; // Start at perigee

        // Calculate position and velocity vectors
        std::vector<double> position(3);
        std::vector<double> velocity(3);
        calculate_orbit_state_vectors(semi_major_axis, eccentricity, inclination, argument_of_perigee, right_ascension, true_anomaly, position, velocity);

        // Generate quaternion for initial attitude
        std::vector<double> quaternion = {0.0, 0.0, 0.0, 1.0};

        return {
            0.0, 0.0, 0.0,  // Angular velocities
            quaternion[0], quaternion[1], quaternion[2], quaternion[3],  // Quaternion
            position[0], position[1], position[2],  // Position
            velocity[0], velocity[1], velocity[2]  // Velocity
        };
    }

    void SatelliteModel::calculate_orbit_state_vectors(double a, double e, double i, double omega, double Omega, double nu,
                                       std::vector<double>& position, std::vector<double>& velocity) {
        double p = a * (1 - e * e);
        double r = p / (1 + e * cos(nu));

        // Position in orbital plane
        double x = r * cos(nu);
        double y = r * sin(nu);

        // Velocity in orbital plane
        double h = sqrt(G * M_earth * p);
        double v_x = -(h / p) * sin(nu);
        double v_y = (h / p) * (e + cos(nu));

        // Rotation matrices
        std::vector<std::vector<double>> R_z = {{cos(Omega), -sin(Omega), 0},
                                                {sin(Omega), cos(Omega), 0},
                                                {0, 0, 1}};
        std::vector<std::vector<double>> R_x = {{1, 0, 0},
                                                {0, cos(i), -sin(i)},
                                                {0, sin(i), cos(i)}};
        std::vector<std::vector<double>> R_z2 = {{cos(omega), -sin(omega), 0},
                                                 {sin(omega), cos(omega), 0},
                                                 {0, 0, 1}};

        // Combine rotation matrices
        auto R = matrix_multiply(R_z, matrix_multiply(R_x, R_z2));

        // Transform position and velocity to Earth-centered inertial frame
        std::vector<double> pos_orbital = {x, y, 0};
        std::vector<double> vel_orbital = {v_x, v_y, 0};

        position = matrix_vector_multiply(R, pos_orbital);
        velocity = matrix_vector_multiply(R, vel_orbital);
    }


    std::vector<double> calculate_moon_position(double t) {
        double theta = n_moon * t;
        double x_moon = a_moon * cos(theta);
        double y_moon = a_moon * sin(theta);
        double z_moon = 0.0; // Assuming Moon's orbit is in the x-y plane
        return {x_moon, y_moon, z_moon};
    }

    std::vector<std::vector<double>> calculate_planet_positions(double t) {
        // Simple circular orbit approximation
        std::vector<std::vector<double>> positions;
        std::vector<double> orbital_periods = {
            0.241 * 365.25 * 24 * 3600, // Mercury
            0.615 * 365.25 * 24 * 3600, // Venus
            1.0 * 365.25 * 24 * 3600,   // Earth
            1.881 * 365.25 * 24 * 3600, // Mars
            11.86 * 365.25 * 24 * 3600, // Jupiter
            29.46 * 365.25 * 24 * 3600, // Saturn
            84.01 * 365.25 * 24 * 3600, // Uranus
            164.8 * 365.25 * 24 * 3600  // Neptune
        };
        std::vector<double> distances = {R_mercury, R_venus, R_earth, R_mars, R_jupiter, R_saturn, R_uranus, R_neptune};

        for (size_t i = 0; i < orbital_periods.size(); ++i) {
            double angle = 2 * M_PI * std::fmod(t, orbital_periods[i]) / orbital_periods[i];
            double x = distances[i] * std::cos(angle);
            double y = distances[i] * std::sin(angle);
            positions.push_back({x, y, 0}); // Assuming orbits are in the x-y plane
        }

        positions.push_back({0, 0, 0}); // Sun at the center
        return positions;
    }

    std::vector<double> SatelliteModel::calculate_lunar_perturbation(double a, double e, double i, double Omega, double omega, double nu, double t) {
        std::vector<double> satellite_pos = keplerian_to_cartesian(a, e, i, Omega, omega, nu);
        std::vector<double> moon_pos = calculate_moon_position(t);
        std::vector<double> r_sat = {satellite_pos[0], satellite_pos[1], satellite_pos[2]};
        std::vector<double> r_moon = {moon_pos[0], moon_pos[1], moon_pos[2]};
        std::vector<double> r_rel = subtract_vectors(r_moon, r_sat);

        double r_rel_mag = magnitude(r_rel);
        double r_moon_mag = magnitude(r_moon);

        std::vector<double> accel = subtract_vectors(
            scale_vector(r_rel, G * M_moon / (r_rel_mag * r_rel_mag * r_rel_mag)),
            scale_vector(r_moon, G * M_moon / (r_moon_mag * r_moon_mag * r_moon_mag))
            );

        return cartesian_to_keplerian_rates(a, e, i, Omega, omega, nu, accel);
    }

    std::vector<double> SatelliteModel::calculate_solar_planetary_perturbations(double a, double e, double i, double Omega, double omega, double nu, double t) {
        std::vector<double> satellite_pos = keplerian_to_cartesian(a, e, i, Omega, omega, nu);
        std::vector<std::vector<double>> body_positions = calculate_planet_positions(t);
        std::vector<double> masses = {M_mercury, M_venus, M_mars, M_jupiter, M_saturn, M_uranus, M_neptune, M_sun};

        std::vector<double> total_accel = {0, 0, 0};

        for (size_t j = 0; j < body_positions.size(); ++j) {
            std::vector<double> r_body = body_positions[j];
            std::vector<double> r_rel = subtract_vectors(r_body, satellite_pos);

            double r_rel_mag = magnitude(r_rel);
            double r_body_mag = magnitude(r_body);

            std::vector<double> accel = subtract_vectors(
                scale_vector(r_rel, G * masses[j] / (r_rel_mag * r_rel_mag * r_rel_mag)),
                scale_vector(r_body, G * masses[j] / (r_body_mag * r_body_mag * r_body_mag))
                );

            total_accel = add_vectors(total_accel, accel);
        }

        return cartesian_to_keplerian_rates(a, e, i, Omega, omega, nu, total_accel);
    }

    std::vector<double> SatelliteModel::cartesian_to_keplerian_rates(double a, double e, double i, double Omega, double omega, double nu, const std::vector<double>& accel) {
        double p = a * (1 - e * e);
        double r = p / (1 + e * std::cos(nu));
        double h = std::sqrt(G * M_earth * p);

        double sin_nu = std::sin(nu);
        double cos_nu = std::cos(nu);

        std::vector<double> e_vec = {e * cos_nu, e * sin_nu, 0};
        std::vector<double> h_vec = {0, 0, h};

        std::vector<double> r_vec = scale_vector(e_vec, r / e);
        std::vector<double> v_vec = cross_product(h_vec, e_vec);
        v_vec = scale_vector(v_vec, std::sqrt(G * M_earth / p));

        double a_dot = 2 * a * a / h * (e * sin_nu * accel[0] + p / r * accel[1]);
        double e_dot = 1 / h * ((p * cos_nu - 2 * r * e) * accel[0] - (p + r) * sin_nu * accel[1]);
        double i_dot = r * cos(omega + nu) / h * accel[2];
        double Omega_dot = r * sin(omega + nu) / (h * sin(i)) * accel[2];
        double omega_dot = 1 / (h * e) * (-p * cos_nu * accel[0] + (p + r) * sin_nu * accel[1]) - r * sin(omega + nu) * cos(i) / (h * sin(i)) * accel[2];
        double nu_dot = h / (r * r) + 1 / (h * e) * (p * cos_nu * accel[0] - (p + r) * sin_nu * accel[1]);

        return {a_dot, e_dot, i_dot, Omega_dot, omega_dot, nu_dot};
    }

    std::vector<double> SatelliteModel::calculate_magnetic_field(double a, double e, double i, double Omega, double omega, double nu, double t) {
        std::vector<double> pos = keplerian_to_cartesian(a, e, i, Omega, omega, nu);
        double x = pos[0], y = pos[1], z = pos[2];
        double r = std::sqrt(x*x + y*y + z*z);
        double theta = std::acos(z / r);
        double phi = std::atan2(y, x);

        double B_r = 0, B_theta = 0, B_phi = 0;

        for (int n = 1; n <= 13; ++n) {
            double r_ratio = a / r;
            double r_ratio_n_plus_1 = std::pow(r_ratio, n + 1);

            for (int m = 0; m <= n; ++m) {
                double P_nm = legendre_polynomial(n, m, std::cos(theta));
                double dP_nm = legendre_polynomial_derivative(n, m, std::cos(theta));

                double cos_m_phi = std::cos(m * phi);
                double sin_m_phi = std::sin(m * phi);

                B_r += (n + 1) * r_ratio_n_plus_1 * (g[n][m] * cos_m_phi + h[n][m] * sin_m_phi) * P_nm;
                B_theta -= r_ratio_n_plus_1 * (g[n][m] * cos_m_phi + h[n][m] * sin_m_phi) * dP_nm;
                B_phi -= m * r_ratio_n_plus_1 * (-g[n][m] * sin_m_phi + h[n][m] * cos_m_phi) * P_nm / std::sin(theta);
            }
        }

        double B_x = B_r * std::sin(theta) * std::cos(phi) + B_theta * std::cos(theta) * std::cos(phi) - B_phi * std::sin(phi);
        double B_y = B_r * std::sin(theta) * std::sin(phi) + B_theta * std::cos(theta) * std::sin(phi) + B_phi * std::cos(phi);
        double B_z = B_r * std::cos(theta) - B_theta * std::sin(theta);

        return {B_x, B_y, B_z};
    }

    std::vector<double> calculate_magnetic_torque(const std::vector<double>& magnetic_moment, const std::vector<double>& magnetic_field) {
        const double epsilon = 1e-10;
        double B_magnitude = std::sqrt(magnetic_field[0]*magnetic_field[0] + magnetic_field[1]*magnetic_field[1] + magnetic_field[2]*magnetic_field[2]);
        if (B_magnitude < epsilon) {
            return {0.0, 0.0, 0.0};
        }
        return {
            magnetic_moment[1] * magnetic_field[2] - magnetic_moment[2] * magnetic_field[1],
            magnetic_moment[2] * magnetic_field[0] - magnetic_moment[0] * magnetic_field[2],
            magnetic_moment[0] * magnetic_field[1] - magnetic_moment[1] * magnetic_field[0]
        };
    }

    std::vector<double> calculate_solar_perturbations(double x, double y, double z, double t) {
        std::vector<std::vector<double>> body_positions = calculate_planet_positions(t);
        std::vector<double> masses = {M_mercury, M_venus, M_mars, M_jupiter, M_saturn, M_uranus, M_neptune, M_sun};
    }

    double SatelliteModel::legendre_polynomial(int n, int m, double x) {
        if (m == 0) {
            if (n == 0) return 1;
            if (n == 1) return x;
            double p0 = 1, p1 = x;
            for (int i = 2; i <= n; ++i) {
                double p2 = ((2 * i - 1) * x * p1 - (i - 1) * p0) / i;
                p0 = p1;
                p1 = p2;
            }
            return p1;
        } else {
            double pmm = 1;
            if (m > 0) {
                double somx2 = std::sqrt((1 - x) * (1 + x));
                double fact = 1;
                for (int i = 1; i <= m; ++i) {
                    pmm *= -fact * somx2;
                    fact += 2;
                }
            }
            if (n == m) return pmm;
            double pmmp1 = x * (2 * m + 1) * pmm;
            if (n == m + 1) return pmmp1;
            double pnm = 0;
            for (int i = m + 2; i <= n; ++i) {
                pnm = ((2 * i - 1) * x * pmmp1 - (i + m - 1) * pmm) / (i - m);
                pmm = pmmp1;
                pmmp1 = pnm;
            }
            return pnm;
        }
    }

    double SatelliteModel::legendre_polynomial_derivative(int n, int m, double x) {
        if (n == m) return n * x * legendre_polynomial(n, n, x) / (x * x - 1);
        return (n * x * legendre_polynomial(n, m, x) - (n + m) * legendre_polynomial(n - 1, m, x)) / (x * x - 1);
    }

    std::vector<double> SatelliteModel::calculate_relativistic_corrections(double a, double e, double i, double Omega, double omega, double nu) {
        std::vector<double> state_cartesian = keplerian_to_cartesian(a, e, i, Omega, omega, nu);
        std::vector<double> position = {state_cartesian[0], state_cartesian[1], state_cartesian[2]};
        std::vector<double> velocity = {state_cartesian[3], state_cartesian[4], state_cartesian[5]};
        double c = 299792458.0; // Speed of light in m/s
        double r = std::sqrt(position[0]*position[0] + position[1]*position[1] + position[2]*position[2]);
        double v_squared = velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2];
        double scalar_correction = G * M_earth / (c * c * r * r * r);

        std::vector<double> corrections(3);
        for (int i = 0; i < 3; ++i) {
            corrections[i] = scalar_correction * (
                                 4 * G * M_earth / r * position[i] -
                                 v_squared * position[i] +
                                 4 * (position[0]*velocity[0] + position[1]*velocity[1] + position[2]*velocity[2]) * velocity[i]
                                 );
        }
        // std::cout << "Relativistic corrections: " << corrections[0] << ", " << corrections[1] << ", " << corrections[2] << std::endl;
        return cartesian_to_keplerian_rates(a, e, i, Omega, omega, nu, corrections);
    }

    std::vector<std::vector<double>> SatelliteModel::matrix_multiply(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) {
        int m = A.size();
        int n = B[0].size();
        int p = A[0].size();

        std::vector<std::vector<double>> C(m, std::vector<double>(n, 0.0));

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < p; ++k) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }

        return C;
    }

    std::vector<double> SatelliteModel::matrix_vector_multiply(const std::vector<std::vector<double>>& A, const std::vector<double>& v) {
        int m = A.size();
        int n = v.size();

        std::vector<double> result(m, 0.0);

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                result[i] += A[i][j] * v[j];
            }
        }

        return result;
    }
