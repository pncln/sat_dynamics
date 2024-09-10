// Implement higher-order gravitational harmonics: Extend beyond J2 to include J3, J4, and higher terms for more accurate gravitational modeling.

// Enhance the atmospheric drag model: Add a sophisticated atmospheric drag calculation, considering factors like solar activity, satellite geometry, and varying atmospheric density.

// Improve the solar radiation pressure model: Implement a detailed solar radiation pressure model that accounts for satellite geometry, reflectivity, and shadowing effects.

// Add magnetic field interactions: Incorporate a high-fidelity Earth magnetic field model and simulate its effects on the satellite, especially if it has magnetic torquers.

// Include relativistic effects: For very precise orbit determination, add relativistic corrections to the equations of motion.

// Enhance third-body perturbations: Improve the lunar and planetary perturbation models by using more accurate ephemeris data and including additional celestial bodies.

// Add thermal effects: Model how temperature changes affect the satellite's structure and components, potentially causing slight changes in its moments of inertia or introducing thermal stresses.

// Implement attitude control systems: Add models for reaction wheels, magnetorquers, or thrusters, and implement control laws for attitude maintenance or maneuvering.

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <array>

// Satellite parameters
const double Ix = 10.0; // Moment of inertia around x-axis
const double Iy = 15.0; // Moment of inertia around y-axis
const double Iz = 20.0; // Moment of inertia around z-axis

// Earth parameters
const double G = 6.67430e-11; // Gravitational constant (m^3 kg^-1 s^-2)
const double M_earth = 5.97e24; // Mass of Earth (kg)
const double R_earth = 6371000; // Radius of Earth (m)
const double J2 = 1.08263e-3; // Earth's J2 coefficient

// Lunar parameters
const double M_moon = 7.34767309e22; // Mass of the Moon (kg)
const double a_moon = 384400000.0; // Semi-major axis of Moon's orbit (m)
const double n_moon = 2.6617e-6; // Mean motion of the Moon (rad/s)

// Solar System parameters
const double M_sun = 1.989e30;    // Mass of the Sun (kg)
const double M_mercury = 3.301e23; // Mass of Mercury (kg)
const double M_venus = 4.867e24;   // Mass of Venus (kg)
const double M_mars = 6.417e23;    // Mass of Mars (kg)
const double M_jupiter = 1.898e27; // Mass of Jupiter (kg)
const double M_saturn = 5.683e26;  // Mass of Saturn (kg)
const double M_uranus = 8.681e25;  // Mass of Uranus (kg)
const double M_neptune = 1.024e26; // Mass of Neptune (kg)

// Average distances from the Sun (m)
const double R_mercury = 5.79e10;
const double R_venus = 1.082e11;
const double R_mars = 2.279e11;
const double R_jupiter = 7.786e11;
const double R_saturn = 1.434e12;
const double R_uranus = 2.871e12;
const double R_neptune = 4.495e12;

// IGRF-13 coefficients (example values, you should use the full set)
// IGRF-13 coefficients (epoch 2020.0)
const std::array<std::array<double, 13>, 13> g = {{
    {0.0, -29404.8, -1450.9, -2500.0, 2982.0, 1676.7, 1363.3, -2381.2, 1236.2, 525.8, 903.0, 809.5, -217.6},
    {0.0, -2499.6, 2982.0, 1677.0, -2991.6, -734.6, 1255.9, 271.5, -231.1, -165.8, 357.4, 68.2, 68.6},
    {0.0, 0.0, -2991.6, -734.8, 1685.9, -575.4, 1251.7, 120.9, -335.6, -86.3, 248.0, -228.0, 92.8},
    {0.0, 0.0, 0.0, 957.1, 245.0, -538.4, 359.6, -157.5, 186.6, -92.9, -109.0, -23.9, 72.9},
    {0.0, 0.0, 0.0, 0.0, 290.2, -227.0, -297.5, 95.2, -30.1, -139.6, -118.2, 78.0, 45.8},
    {0.0, 0.0, 0.0, 0.0, 0.0, 376.0, 47.7, 275.9, -178.3, -51.4, -6.8, -16.9, -2.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 67.6, -25.9, -1.2, 26.4, 17.1, -50.8, 14.7},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -26.3, 4.6, 24.2, 8.9, 10.1, -14.1},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.3, -6.2, -18.3, 7.0, 9.4},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -10.4, -3.4, 2.4, -0.4},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, -6.5, 5.7},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.6, -0.5},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1}
}};

const std::array<std::array<double, 13>, 13> h = {{
    {0.0, 0.0, 4652.5, 0.0, -2991.6, -734.6, 1255.9, 271.5, -231.1, -165.8, 357.4, 68.2, 68.6},
    {0.0, 0.0, 0.0, -734.8, 1685.9, -575.4, 1251.7, 120.9, -335.6, -86.3, 248.0, -228.0, 92.8},
    {0.0, 0.0, 0.0, 0.0, 245.0, -538.4, 359.6, -157.5, 186.6, -92.9, -109.0, -23.9, 72.9},
    {0.0, 0.0, 0.0, 0.0, 0.0, -227.0, -297.5, 95.2, -30.1, -139.6, -118.2, 78.0, 45.8},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 47.7, 275.9, -178.3, -51.4, -6.8, -16.9, -2.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -25.9, -1.2, 26.4, 17.1, -50.8, 14.7},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.6, 24.2, 8.9, 10.1, -14.1},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -6.2, -18.3, 7.0, 9.4},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.4, 2.4, -0.4},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -6.5, 5.7},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
}};


const double a = 6371.2; // Earth's mean radius in km

// State space model for satellite rotational dynamics
class SatelliteModel {
private:
    std::vector<double> state; // [wx, wy, wz, qx, qy, qz, q0, x, y, z, vx, vy, vz]
    double mass; // Mass of the satellite (kg)
    std::vector<double> magnetic_moment;

public:
    SatelliteModel(const std::vector<double>& initial_state, double satellite_mass, const std::vector<double>& initial_magnetic_moment) : state(initial_state), mass(satellite_mass), magnetic_moment(initial_magnetic_moment) {}

    void update(double dt, const std::vector<double>& torque, const std::vector<double>& external_disturbance, double t) {
        // Rotational dynamics equations in state space form
        double wx = state[0];
        double wy = state[1];
        double wz = state[2];

        // Position of the satellite
        double x = state[7];
        double y = state[8];
        double z = state[9];

        // Velocity of the satellite
        double vx = state[10];
        double vy = state[11];
        double vz = state[12];

        // J2 perturbation
        std::vector<double> j2_perturbation = calculate_j2_perturbation(x, y, z);

        // gravitational force torque
        std::vector<double> grav_force = calculate_gravitational_force(x, y, z);
        double F_grav = std::sqrt(grav_force[0]*grav_force[0] + grav_force[1]*grav_force[1] + grav_force[2]*grav_force[2]);
        std::vector<double> grav_torque = calculate_gravitational_torque(x, y, z, F_grav);

        // lunar perturbation
        std::vector<double> lunar_pert = calculate_lunar_perturbation(x, y, z, t);

        std::vector<double> gyro_torque = calculate_gyroscopic_torque(wx, wy, wz);

        double wx_dot = (torque[0] + external_disturbance[0] + grav_torque[0] + gyro_torque[0] + (Iy - Iz) * wy * wz) / Ix;
        double wy_dot = (torque[1] + external_disturbance[1] + grav_torque[1] + gyro_torque[1] + (Iz - Ix) * wz * wx) / Iy;
        double wz_dot = (torque[2] + external_disturbance[2] + grav_torque[2] + gyro_torque[2] + (Ix - Iy) * wx * wy) / Iz;

        std::vector<double> k1 = {wx_dot, wy_dot, wz_dot};
        std::vector<double> k2 = calculate_derivatives(wx + 0.5 * dt * k1[0], wy + 0.5 * dt * k1[1], wz + 0.5 * dt * k1[2], torque, external_disturbance, x, y, z);
        std::vector<double> k3 = calculate_derivatives(wx + 0.5 * dt * k2[0], wy + 0.5 * dt * k2[1], wz + 0.5 * dt * k2[2], torque, external_disturbance, x, y, z);
        std::vector<double> k4 = calculate_derivatives(wx + dt * k3[0], wy + dt * k3[1], wz + dt * k3[2], torque, external_disturbance, x, y, z);

        state[0] += (dt / 6.0) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
        state[1] += (dt / 6.0) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
        state[2] += (dt / 6.0) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);

        // Quaternion kinematics
        double qx = state[3];
        double qy = state[4];
        double qz = state[5];
        double q0 = state[6];

        // Quaternion derivative
        double qx_dot = 0.5 * (q0 * wx - qz * wy + qy * wz);
        double qy_dot = 0.5 * (qz * wx + q0 * wy - qx * wz);
        double qz_dot = 0.5 * (-qy * wx + qx * wy + q0 * wz);
        double q0_dot = -0.5 * (qx * wx + qy * wy + qz * wz);

        std::vector<double> q_k1 = {qx_dot, qy_dot, qz_dot, q0_dot};
        std::vector<double> q_k2 = calculate_quaternion_derivatives(qx + 0.5 * dt * q_k1[0], qy + 0.5 * dt * q_k1[1], qz + 0.5 * dt * q_k1[2], q0 + 0.5 * dt * q_k1[3], wx, wy, wz);
        std::vector<double> q_k3 = calculate_quaternion_derivatives(qx + 0.5 * dt * q_k2[0], qy + 0.5 * dt * q_k2[1], qz + 0.5 * dt * q_k2[2], q0 + 0.5 * dt * q_k2[3], wx, wy, wz);
        std::vector<double> q_k4 = calculate_quaternion_derivatives(qx + dt * q_k3[0], qy + dt * q_k3[1], qz + dt * q_k3[2], q0 + dt * q_k3[3], wx, wy, wz);

        state[3] += (dt / 6.0) * (q_k1[0] + 2 * q_k2[0] + 2 * q_k3[0] + q_k4[0]);
        state[4] += (dt / 6.0) * (q_k1[1] + 2 * q_k2[1] + 2 * q_k3[1] + q_k4[1]);
        state[5] += (dt / 6.0) * (q_k1[2] + 2 * q_k2[2] + 2 * q_k3[2] + q_k4[2]);
        state[6] += (dt / 6.0) * (q_k1[3] + 2 * q_k2[3] + 2 * q_k3[3] + q_k4[3]);

        double norm = sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5] + state[6]*state[6]);
        state[3] /= norm;
        state[4] /= norm;
        state[5] /= norm;
        state[6] /= norm;

        // Calculate magnetic field at the satellite's position
        std::vector<double> magnetic_field = calculate_magnetic_field(x, y, z, t);

        // Calculate magnetic torque
        std::vector<double> magnetic_torque = calculate_magnetic_torque(magnetic_moment, magnetic_field);

        // Add magnetic torque to the rotational dynamics
        wx_dot += magnetic_torque[0] / Ix;
        wy_dot += magnetic_torque[1] / Iy;
        wz_dot += magnetic_torque[2] / Iz;

        std::vector<double> solar_planetary_pert = calculate_solar_planetary_perturbations(x, y, z, t);

        std::vector<double> r_k1 = {vx, vy, vz};
        std::vector<double> v_k1 = {(grav_force[0] + j2_perturbation[0] + lunar_pert[0] + solar_planetary_pert[0]) / mass, 
                                    (grav_force[1] + j2_perturbation[1] + lunar_pert[1] + solar_planetary_pert[1]) / mass, 
                                    (grav_force[2] + j2_perturbation[2] + lunar_pert[2] + solar_planetary_pert[2]) / mass};

        std::vector<double> r_k2 = {vx + 0.5 * dt * v_k1[0], vy + 0.5 * dt * v_k1[1], vz + 0.5 * dt * v_k1[2]};
        std::vector<double> grav_force_k2 = calculate_gravitational_force(x + 0.5 * dt * r_k1[0], y + 0.5 * dt * r_k1[1], z + 0.5 * dt * r_k1[2]);
        std::vector<double> j2_perturbation_k2 = calculate_j2_perturbation(x + 0.5 * dt * r_k1[0], y + 0.5 * dt * r_k1[1], z + 0.5 * dt * r_k1[2]);
        std::vector<double> solar_planetary_pert_k2 = calculate_solar_planetary_perturbations(x + 0.5 * dt * r_k1[0], y + 0.5 * dt * r_k1[1], z + 0.5 * dt * r_k1[2], t + 0.5 * dt);
        std::vector<double> lunar_pert_k2 = calculate_lunar_perturbation(x + 0.5 * dt * r_k1[0], y + 0.5 * dt * r_k1[1], z + 0.5 * dt * r_k1[2], t + 0.5 * dt);
        std::vector<double> v_k2 = {(grav_force_k2[0] + j2_perturbation_k2[0] + lunar_pert_k2[0] + solar_planetary_pert_k2[0]) / mass, 
                                    (grav_force_k2[1] + j2_perturbation_k2[1] + lunar_pert_k2[1] + solar_planetary_pert_k2[1]) / mass, 
                                    (grav_force_k2[2] + j2_perturbation_k2[2] + lunar_pert_k2[2] + solar_planetary_pert_k2[2]) / mass};

        std::vector<double> r_k3 = {vx + 0.5 * dt * v_k2[0], vy + 0.5 * dt * v_k2[1], vz + 0.5 * dt * v_k2[2]};
        std::vector<double> grav_force_k3 = calculate_gravitational_force(x + 0.5 * dt * r_k2[0], y + 0.5 * dt * r_k2[1], z + 0.5 * dt * r_k2[2]);
        std::vector<double> j2_perturbation_k3 = calculate_j2_perturbation(x + 0.5 * dt * r_k2[0], y + 0.5 * dt * r_k2[1], z + 0.5 * dt * r_k2[2]);
        std::vector<double> solar_planetary_pert_k3 = calculate_solar_planetary_perturbations(x + 0.5 * dt * r_k2[0], y + 0.5 * dt * r_k2[1], z + 0.5 * dt * r_k2[2], t + 0.5 * dt);
        std::vector<double> lunar_pert_k3 = calculate_lunar_perturbation(x + 0.5 * dt * r_k2[0], y + 0.5 * dt * r_k2[1], z + 0.5 * dt * r_k2[2], t + 0.5 * dt);
        std::vector<double> v_k3 = {(grav_force_k3[0] + j2_perturbation_k3[0] + lunar_pert_k3[0] + solar_planetary_pert_k3[0]) / mass, 
                                    (grav_force_k3[1] + j2_perturbation_k3[1] + lunar_pert_k3[1] + solar_planetary_pert_k3[1]) / mass, 
                                    (grav_force_k3[2] + j2_perturbation_k3[2] + lunar_pert_k3[2] + solar_planetary_pert_k3[2]) / mass};

        std::vector<double> r_k4 = {vx + dt * v_k3[0], vy + dt * v_k3[1], vz + dt * v_k3[2]};
        std::vector<double> grav_force_k4 = calculate_gravitational_force(x + dt * r_k3[0], y + dt * r_k3[1], z + dt * r_k3[2]);
        std::vector<double> j2_perturbation_k4 = calculate_j2_perturbation(x + dt * r_k3[0], y + dt * r_k3[1], z + dt * r_k3[2]);
        std::vector<double> solar_planetary_pert_k4 = calculate_solar_planetary_perturbations(x + dt * r_k3[0], y + dt * r_k3[1], z + dt * r_k3[2], t + dt);
        std::vector<double> lunar_pert_k4 = calculate_lunar_perturbation(x + dt * r_k3[0], y + dt * r_k3[1], z + dt * r_k3[2], t + dt);
        std::vector<double> v_k4 = {(grav_force_k4[0] + j2_perturbation_k4[0] + lunar_pert_k4[0] + solar_planetary_pert_k4[0]) / mass, 
                                    (grav_force_k4[1] + j2_perturbation_k4[1] + lunar_pert_k4[1] + solar_planetary_pert_k4[1]) / mass, 
                                    (grav_force_k4[2] + j2_perturbation_k4[2] + lunar_pert_k4[2] + solar_planetary_pert_k4[2]) / mass};
        

        state[7] += (dt / 6.0) * (r_k1[0] + 2 * r_k2[0] + 2 * r_k3[0] + r_k4[0]);
        state[8] += (dt / 6.0) * (r_k1[1] + 2 * r_k2[1] + 2 * r_k3[1] + r_k4[1]);
        state[9] += (dt / 6.0) * (r_k1[2] + 2 * r_k2[2] + 2 * r_k3[2] + r_k4[2]);

        state[10] += (dt / 6.0) * (v_k1[0] + 2 * v_k2[0] + 2 * v_k3[0] + v_k4[0]);
        state[11] += (dt / 6.0) * (v_k1[1] + 2 * v_k2[1] + 2 * v_k3[1] + v_k4[1]);
        state[12] += (dt / 6.0) * (v_k1[2] + 2 * v_k2[2] + 2 * v_k3[2] + v_k4[2]);
    }

    const std::vector<double>& get_state() const {
        return state;
    }

private:
    std::vector<double> calculate_derivatives(double wx, double wy, double wz, const std::vector<double>& torque, const std::vector<double>& external_disturbance, double x, double y, double z) {
        std::vector<double> grav_force = calculate_gravitational_force(x, y, z);
        double F_grav = std::sqrt(grav_force[0]*grav_force[0] + grav_force[1]*grav_force[1] + grav_force[2]*grav_force[2]);
        std::vector<double> grav_torque = calculate_gravitational_torque(x, y, z, F_grav);
        std::vector<double> gyro_torque = calculate_gyroscopic_torque(wx, wy, wz);

        double wx_dot = (torque[0] + external_disturbance[0] + grav_torque[0] + gyro_torque[0] + (Iy - Iz) * wy * wz) / Ix;
        double wy_dot = (torque[1] + external_disturbance[1] + grav_torque[1] + gyro_torque[1] + (Iz - Ix) * wz * wx) / Iy;
        double wz_dot = (torque[2] + external_disturbance[2] + grav_torque[2] + gyro_torque[2] + (Ix - Iy) * wx * wy) / Iz;
        return {wx_dot, wy_dot, wz_dot};
    }
    std::vector<double> calculate_quaternion_derivatives(double qx, double qy, double qz, double q0, double wx, double wy, double wz) {
        double qx_dot = 0.5 * (q0 * wx - qz * wy + qy * wz);
        double qy_dot = 0.5 * (qz * wx + q0 * wy - qx * wz);
        double qz_dot = 0.5 * (-qy * wx + qx * wy + q0 * wz);
        double q0_dot = -0.5 * (qx * wx + qy * wy + qz * wz);
        return {qx_dot, qy_dot, qz_dot, q0_dot};
    }
    std::vector<double> calculate_gravitational_torque(double x, double y, double z, double F_grav) {
        double r = std::sqrt(x*x + y*y + z*z);
        double Tx = 3 * F_grav * (y*z/r) * (Iz - Iy) / r;
        double Ty = 3 * F_grav * (z*x/r) * (Ix - Iz) / r;
        double Tz = 3 * F_grav * (x*y/r) * (Iy - Ix) / r;
        return {Tx, Ty, Tz};
    }
    std::vector<double> calculate_gyroscopic_torque(double wx, double wy, double wz) {
        double Tx = (Iy - Iz) * wy * wz;
        double Ty = (Iz - Ix) * wz * wx;
        double Tz = (Ix - Iy) * wx * wy;
        return {Tx, Ty, Tz};
    }

    std::vector<double> calculate_gravitational_force(double x, double y, double z) {
        double r = std::sqrt(x*x + y*y + z*z);
        double F = -G * M_earth * mass / (r * r);
        return {F * x / r, F * y / r, F * z / r};
    }

    std::vector<double> calculate_j2_perturbation(double x, double y, double z) {
        double r = std::sqrt(x*x + y*y + z*z);
        double r_squared = r * r;
        double r_cubed = r * r_squared;
        double r_fifth = r_cubed * r_squared;
        double coeff = -1.5 * J2 * G * M_earth * R_earth * R_earth / r_fifth;
        
        double x_pert = coeff * x * (5 * z * z / r_squared - 1);
        double y_pert = coeff * y * (5 * z * z / r_squared - 1);
        double z_pert = coeff * z * (5 * z * z / r_squared - 3);
        
        return {x_pert, y_pert, z_pert};
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

    std::vector<double> calculate_lunar_perturbation(double x, double y, double z, double t) {
        std::vector<double> moon_pos = calculate_moon_position(t);
        double dx = x - moon_pos[0];
        double dy = y - moon_pos[1];
        double dz = z - moon_pos[2];
        double r_moon = std::sqrt(dx*dx + dy*dy + dz*dz);
        double r_moon_cubed = r_moon * r_moon * r_moon;
        
        double ax = G * M_moon * (moon_pos[0] / (a_moon * a_moon * a_moon) - dx / r_moon_cubed);
        double ay = G * M_moon * (moon_pos[1] / (a_moon * a_moon * a_moon) - dy / r_moon_cubed);
        double az = G * M_moon * (moon_pos[2] / (a_moon * a_moon * a_moon) - dz / r_moon_cubed);
        
        return {ax, ay, az};
    }

    std::vector<double> calculate_solar_planetary_perturbations(double x, double y, double z, double t) {
        std::vector<std::vector<double>> body_positions = calculate_planet_positions(t);
        std::vector<double> masses = {M_mercury, M_venus, M_mars, M_jupiter, M_saturn, M_uranus, M_neptune, M_sun};

        double ax = 0, ay = 0, az = 0;

        for (size_t i = 0; i < body_positions.size(); ++i) {
            double dx = x - body_positions[i][0];
            double dy = y - body_positions[i][1];
            double dz = z - body_positions[i][2];
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            double r_cubed = r * r * r;

            ax += G * masses[i] * (-dx / r_cubed);
            ay += G * masses[i] * (-dy / r_cubed);
            az += G * masses[i] * (-dz / r_cubed);
        }

        return {ax, ay, az};
    }

    std::vector<double> calculate_magnetic_field(double x, double y, double z, double t) {
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

    double legendre_polynomial(int n, int m, double x) {
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

    double legendre_polynomial_derivative(int n, int m, double x) {
        if (n == m) return n * x * legendre_polynomial(n, n, x) / (x * x - 1);
        return (n * x * legendre_polynomial(n, m, x) - (n + m) * legendre_polynomial(n - 1, m, x)) / (x * x - 1);
    }
};

int main() {
    // Initial state: [wx, wy, wz, qx, qy, qz, q0, x, y, z, vx, vy, vz]
    std::vector<double> initial_state = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, R_earth + 500000, 0.0, 0.0, 0.0, std::sqrt(G * M_earth / (R_earth + 500000)), 0.0};
    double satellite_mass = 1000.0; // kg
    std::vector<double> initial_magnetic_moment = {1.0, 1.0, 1.0}; // A·m² (example values)

    double dt = 0.001; // Time step
    double t_max = 10; // Total simulation time

    // Impulse thrust/torque
    double impulse_duration = 1.0; // Duration of impulse thrust/torque
    std::vector<double> impulse_torque = {0.1, 0.1, 0.1}; // Torque applied during impulse

    // Prepare data for plotting
    std::vector<double> time_data;
    std::vector<std::vector<double>> state_data(13);

    SatelliteModel satellite(initial_state, satellite_mass, initial_magnetic_moment);

    // Simulation loop
    for (double t = 0; t < t_max; t += dt) {
        std::vector<double> torque = {0.0, 0.0, 0.0};
        std::vector<double> external_disturbance = {0.0, 0.0, 0.0};
        
        if (t < impulse_duration) {
            torque = impulse_torque;
        }

        satellite.update(dt, torque, external_disturbance, t);

        // Store data for plotting
        time_data.push_back(t);
        const auto& state = satellite.get_state();
        for (int i = 0; i < 13; ++i) {
            state_data[i].push_back(state[i]);
        }

        if (std::fmod(t, 1.0) <= dt) {
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "Time: " << std::setw(8) << t << " seconds" << std::endl;
            std::cout << "+--------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+" << std::endl;
            std::cout << "| Param  |     wx     |     wy     |     wz     |     qx     |     qy     |     qz     |     q0     |     x      |     y      |     z      |     vx     |     vy     |     vz     |" << std::endl;
            std::cout << "+--------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+" << std::endl;
            std::cout << "| Value  | " << std::setw(10) << state[0] << " | " << std::setw(10) << state[1] << " | "
                      << std::setw(10) << state[2] << " | " << std::setw(10) << state[3] << " | "
                      << std::setw(10) << state[4] << " | " << std::setw(10) << state[5] << " | "
                      << std::setw(10) << state[6] << " | " << std::setw(10) << state[7] << " | "
                      << std::setw(10) << state[8] << " | " << std::setw(10) << state[9] << " | "
                      << std::setw(10) << state[10] << " | " << std::setw(10) << state[11] << " | "
                      << std::setw(10) << state[12] << " |" << std::endl;
            std::cout << "+--------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+------------+" << std::endl;
        }
    }

    std::ofstream outfile("satellite_data.csv");
    outfile << "Time,wx,wy,wz,qx,qy,qz,q0,x,y,z,vx,vy,vz\n";

    for (size_t i = 0; i < time_data.size(); ++i) {
        outfile << time_data[i];
        for (const auto& param : state_data) {
            outfile << "," << param[i];
        }
        outfile << "\n";
    }
    outfile.close();

    std::cout << "Data has been written to 'satellite_data.csv'." << std::endl;
    std::cout << "Started plotting data..." << std::endl;

    int result = std::system("./plot");

    if (result == 0) {
        std::cout << "External process completed successfully." << std::endl;
    } else {
        std::cout << "External process failed with exit code: " << result << std::endl;
    }

    return 0;
}
