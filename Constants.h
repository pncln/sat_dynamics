#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <array>
#include <cmath>

// Satellite parameters
const double Ix = 10.0; // Moment of inertia around x-axis
const double Iy = 15.0; // Moment of inertia around y-axis
const double Iz = 20.0; // Moment of inertia around z-axis

// Earth parameters
const double G = 6.67430e-11; // Gravitational constant (m^3 kg^-1 s^-2)
const double M_earth = 5.97e24; // Mass of Earth (kg)
const double R_earth = 6371000; // Radius of Earth (m)

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

#endif // CONSTANTS_H
