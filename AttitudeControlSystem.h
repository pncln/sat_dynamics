#ifndef ATTITUDE_CONTROL_SYSTEM_H
#define ATTITUDE_CONTROL_SYSTEM_H

#include <vector>
class AttitudeControlSystem {
private:
    std::vector<double> reaction_wheel_momentum;
    std::vector<double> magnetorquer_dipole;
    std::vector<double> thruster_force;
    std::vector<double> Kp;
    std::vector<double> Kd;

public:
    AttitudeControlSystem();

    std::vector<double> calculate_control_torque(const std::vector<double>& current_state, const std::vector<double>& desired_state, const std::vector<double>& magnetic_field);
    void allocate_control_torque(const std::vector<double>& control_torque, const std::vector<double>& magnetic_field);
    void calculate_magnetorquer_dipole(const std::vector<double>& desired_torque, const std::vector<double>& magnetic_field);

    std::vector<double> get_reaction_wheel_momentum() const;
    std::vector<double> get_magnetorquer_dipole() const;
    std::vector<double> get_thruster_force() const;
};

#endif // ATTITUDE_CONTROL_SYSTEM_H
