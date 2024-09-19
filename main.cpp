#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include "SatelliteModel.h"
#include "AttitudeControlSystem.h"


int main() {
    // Initial state: [wx, wy, wz, qx, qy, qz, q0, a, e, i, Omega, omega, nu]
    std::vector<double> initial_state = {
        0.0, 0.0, 0.0,                 // Angular velocities
        0.0, 0.0, 0.0, 1.0,            // Quaternion (identity rotation)
        R_earth + 500000.0,            // Semi-major axis (500 km altitude)
        0.001,                         // Eccentricity (near-circular orbit)
        0.9553,                        // Inclination (54.74 degrees)
        0.0,                           // Right Ascension of Ascending Node
        0.0,                           // Argument of Perigee
        0.0                            // True Anomaly
    };
    double satellite_mass = 1000.0; // kg
    std::vector<double> initial_magnetic_moment = {1.0, 1.0, 1.0}; // A·m² (example values)

    SatelliteModel satellite(initial_state, satellite_mass, initial_magnetic_moment);
    AttitudeControlSystem acs;

    double dt = 0.05; // Time step
    double t_max = 60*60; // Total simulation time ()

    // Prepare data for plotting
    std::vector<double> time_data;
    std::vector<std::vector<double>> state_data(13);

    std::vector<double> torque = {0.1, 0.1, 0.1};
    std::vector<double> external_disturbance = {0.2, 0.2, 0.2};

    // Simulation loop
    for (double t = 0; t < t_max; t += dt) {
        const auto& current_state = satellite.get_state();
        std::vector<double> desired_state = current_state; // For now, just maintain current state
        std::vector<double> magnetic_field = satellite.calculate_magnetic_field(
            current_state[7], current_state[8], current_state[9],
            current_state[10], current_state[11], current_state[12], t
        );

        std::vector<double> control_torque = acs.calculate_control_torque(current_state, desired_state, magnetic_field);

        satellite.update(dt, torque, external_disturbance, t); // Assuming no external disturbance for simplicity

        // Store data for plotting
        time_data.push_back(t);
        const auto& state = satellite.get_state();
        for (int i = 0; i < 13; ++i) {
            state_data[i].push_back(state[i]);
        }

        // Orbit-related checks
        double orbital_period = 2 * M_PI * std::sqrt(std::pow(state[7], 3) / (G * M_earth));
        if (std::fmod(t, orbital_period) < dt) {
            std::cout << "Completed one orbit at time: " << t << " seconds" << std::endl;
        }

        // Print status every minute of simulation time
        if (std::fmod(t, 60) < dt) {
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "Time: " << std::setw(8) << t / 3600.0 << " hours" << std::endl;
            std::cout << "Semi-major axis: " << state[7] << " m, Eccentricity: " << state[8]
                      << ", Inclination: " << state[9] * 180.0 / M_PI << " deg" << std::endl;
            std::cout << "RAAN: " << state[10] * 180.0 / M_PI << " deg, Arg of Perigee: "
                      << state[11] * 180.0 / M_PI << " deg, True Anomaly: " << state[12] * 180.0 / M_PI << " deg" << std::endl;
            std::cout << std::endl;
        }
    }

    // Output data to CSV file
    std::ofstream outfile("satellite_data.csv");
    outfile << "Time,wx,wy,wz,qx,qy,qz,q0,a,e,i,Omega,omega,nu\n";

    for (size_t i = 0; i < time_data.size(); ++i) {
        outfile << time_data[i];
        for (const auto& param : state_data) {
            outfile << "," << param[i];
        }
        outfile << "\n";
    }
    outfile.close();

    std::cout << "Data has been written to 'satellite_data.csv'." << std::endl;

    return 0;
}
