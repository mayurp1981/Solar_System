# Solar_System

This MATLAB program simulates the Solar System by calculating the initial state vectors of planets on a specified Julian date and visualizes and animates their motion in space. The program accounts for gravitational forces from the Sun and other planets, addressing the n-body problem in Astrodynamics. 

The `UTC2JD` function converts the Gregorian date to Julian date for calculating the initial Keplerian elements of planets. These elements are then processed in the `PlanetStateVec` function to obtain state vectors from Keplerian elements, utilizing the `Kep2States` function for the conversion. The `rates` function calculates accelerations and velocities to propagate the orbit of planets.

RK4 numerical integration is employed in this study due to its modularity and ease of implementation.

**Note:** Keep all files in the same folder for seamless integration.
