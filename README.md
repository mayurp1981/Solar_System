# Solar_System

This MATLAB program simulates the Solar System by calculating the initial state vectors of planets on a specified Julian date and visualizes and animates their motion in space. The program accounts for gravitational forces from the Sun and other planets, addressing the n-body problem in Astrodynamics. 

The `UTC2JD` function converts the Gregorian date to the Julian date for calculating the initial Keplerian elements of planets. These elements are then processed in the `PlanetStateVec` function to obtain state vectors from Keplerian elements, utilizing the `Kep2States` function for the conversion. The `rates` function calculates accelerations and velocities to propagate the orbit of planets.

RK4 numerical integration is employed in this study due to its modularity and ease of implementation.

**Note:** Keep all files in the same folder for seamless integration.

Visualization of the inner Solar System, showcasing the dynamic orbits and positions of inner planets around the Sun.
![Untitled design](https://github.com/mayurp1981/Solar_System/assets/75112223/4c3943c9-f8ab-47b0-a2f9-9762a0776ed3)

Comprehensive simulation of the entire Solar System for 5 years.
![Solar_system_B](https://github.com/mayurp1981/Solar_System/assets/75112223/f24a1c87-8e1c-48ad-b9d2-a74a4eed66f1)
