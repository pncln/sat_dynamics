# Satellite Dynamics Simulation & Orbit Propagation

C++ project for simulating complex satellite dynamics in space, accounting for various perturbations and celestial body influences.

## Project Structure

- `main.cpp`: Core simulation code
- `CMakeLists.txt.user`: Qt Creator project configuration
- `.vscode/settings.json`: VS Code settings for C++ development
- `.gitattributes`: Git attributes for the repository

## Detailed Code Overview

### Main Simulation (`main.cpp`)

The `main.cpp` file is the heart of the simulation, containing detailed models for satellite dynamics in space. It includes:

#### Constants and Parameters

1. Satellite Parameters:
   - Moments of inertia (Ix, Iy, Iz) representing the satellite's mass distribution

2. Earth Parameters:
   - Gravitational constant (G)
   - Earth's mass (M_earth)
   - Earth's radius (R_earth)
   - J2 coefficient for Earth's oblateness

3. Lunar Parameters:
   - Moon's mass (M_moon)
   - Semi-major axis of Moon's orbit (a_moon)
   - Mean motion of the Moon (n_moon)

4. Solar System Parameters:
   - Masses of the Sun and Mercury

#### Simulation Components

The simulation aims to model the following effects on satellite dynamics:

1. Earth's Gravitational Field:
   - Basic gravitational attraction
   - J2 perturbation due to Earth's non-spherical shape

2. Lunar Gravitational Effects:
   - Gravitational influence of the Moon on the satellite's orbit

3. Solar and Planetary Influences:
   - Gravitational effects from the Sun and Mercury

#### Planned Enhancements

The code includes detailed comments outlining future improvements:

1. Atmospheric Drag Modeling:
   - Implementation of drag effects for low Earth orbits
   - Calculation of satellite's cross-sectional area and drag coefficient

2. Solar Radiation Pressure:
   - Modeling the effect of solar radiation, crucial for satellites with large solar panels

3. Magnetic Torque:
   - Incorporation of Earth's magnetic field effects, especially for satellites with magnetic torquers

4. Enhanced Gravity Model:
   - Implementation of higher-order gravitational harmonics (J3, J4) for more accurate modeling

5. Thermal Effects:
   - Modeling of temperature-induced changes in satellite structure and moments of inertia

6. Attitude Control Systems:
   - Addition of reaction wheels, magnetorquers, or thrusters models
   - Implementation of control laws for attitude maintenance and maneuvering

## Development Environment

The project is configured for development in both Qt Creator and Visual Studio Code:

### Qt Creator
- Utilizes CMake for project management
- Configuration details in `CMakeLists.txt.user`

### Visual Studio Code
- Configured for C++ development
- `.vscode/settings.json` contains specific file associations for C++ files

## Version Control

Git is used for version control with the following configuration:

- `.gitattributes` file set up for:
  - Automatic text file detection
  - Line ending normalization (LF)

## Getting Started

To run the satellite dynamics simulation:

1. Clone the repository:
```bash
git clone https://github.com/pncln/sat_dynamics.git
```
2. Open the project:
- In Qt Creator: Open the `CMakeLists.txt` file
- In VS Code: Open the project folder
3. Build the project using the IDE's build system
4. Run the `main.cpp` file to start the simulation

## Contributing

Contributions to enhance the simulation are highly encouraged. Potential areas for contribution include:

- Implementing any of the planned enhancements mentioned in `main.cpp`
- Improving the accuracy of existing models
- Adding visualization capabilities for the simulation results
- Optimizing the code for better performance

Please submit pull requests with clear descriptions of the changes and their purposes.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.

## Contact
@pncln
