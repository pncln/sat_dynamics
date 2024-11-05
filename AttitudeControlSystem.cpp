#include <iostream>
#include "AttitudeControlSystem.h"

AttitudeControlSystem::AttitudeControlSystem() :
        reaction_wheel_momentum(3.0, 0.0),
        magnetorquer_dipole(3.0, 0.0),
        thruster_force(3.0, 0.0),
        Kp({0.1, 0.1, 0.05}),
        Kd({0.5, 0.5, 0.25}) {}

std::vector<double> AttitudeControlSystem::calculate_control_torque(const std::vector<double>& current_state, const std::vector<double>& desired_state, const std::vector<double>& magnetic_field) {
        std::vector<double> control_torque(3.0, 0.0);

        for (int i = 0; i < 3; ++i) {
            double error = desired_state[i] - current_state[i];
            double error_rate = -current_state[i]; // Assuming angular velocity starts at index 0
            control_torque[i] = Kp[i] * error + Kd[i] * error_rate;
        }

        // Allocate control torque to actuators
        allocate_control_torque(control_torque, magnetic_field);

        // std::cout << "Control torque: " << control_torque[0] << ", " << control_torque[1] << ", " << control_torque[2] << std::endl;

        return control_torque;
    }

    void AttitudeControlSystem::allocate_control_torque(const std::vector<double>& control_torque, const std::vector<double>& magnetic_field) {
        // Allocate torque to reaction wheels
        for (int i = 0; i < 3; ++i) {
            reaction_wheel_momentum[i] += control_torque[i] * 0.5; // Assume 50% allocation to reaction wheels
        }

        // Allocate torque to magnetorquers
        std::vector<double> magnetorquer_torque(3, 0.0);
        for (int i = 0; i < 3; ++i) {
            magnetorquer_torque[i] = control_torque[i] * 0.3; // Assume 30% allocation to magnetorquers
        }
        calculate_magnetorquer_dipole(magnetorquer_torque, magnetic_field);

        // Allocate remaining torque to thrusters
        for (int i = 0; i < 3; ++i) {
            thruster_force[i] = control_torque[i] * 0.2 / 0.1; // Assume 20% allocation and 0.1m moment arm
        }

        // std::cout << "Reaction wheel momentum: " << reaction_wheel_momentum[0] << ", " << reaction_wheel_momentum[1] << ", " << reaction_wheel_momentum[2] << std::endl;
        // std::cout << "Magnetorquer dipole: " << magnetorquer_dipole[0] << ", " << magnetorquer_dipole[1] << ", " << magnetorquer_dipole[2] << std::endl;
        // std::cout << "Thruster force: " << thruster_force[0] << ", " << thruster_force[1] << ", " << thruster_force[2] << std::endl;
    }

    void AttitudeControlSystem::calculate_magnetorquer_dipole(const std::vector<double>& desired_torque, const std::vector<double>& magnetic_field) {
        // Calculate required magnetic dipole moment
        double B_magnitude = std::sqrt(magnetic_field[0]*magnetic_field[0] + magnetic_field[1]*magnetic_field[1] + magnetic_field[2]*magnetic_field[2]);
        for (int i = 0; i < 3; ++i) {
            magnetorquer_dipole[i] = desired_torque[i] / B_magnitude;
        }
    }

    std::vector<double> AttitudeControlSystem::get_reaction_wheel_momentum() const { return reaction_wheel_momentum; }
    std::vector<double> AttitudeControlSystem::get_magnetorquer_dipole() const { return magnetorquer_dipole; }
    std::vector<double> AttitudeControlSystem::get_thruster_force() const { return thruster_force; }
