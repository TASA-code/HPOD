# Orbit_Determination

This code serves the purpose for orbit determination of a satellite orbiting aroung Earth with the effect of J2 perturbation. 
The code takes classical orbital element (SMA, e, i, RAAN, &omega, &nu) into account to calculate and iterate the final position and velocity vector

- SMA  : Semi-Major axis,
- e    : Eccentricity,
- i    : Inclination,
- RAAN : Right Ascension of Ascending Node,
- &omega: Argument of periapsi,
- &nu : True Anomaly.


For time-integration we will use 4th-order Runge-Kutta explicit scheme. This can be formulated as:

$$
y_{n+1} = y_n + \frac{1}{6}(k_1+2k_2+2k_3+k_4)\Delta t
$$

with 

$$
\begin{align}
    k_1 = f(y_n)\\
    k_2 = f(y_n + \Delta tk_1/2)\\
    k_3 = f(y_n + \Delta tk_2/2)\\
    k_4 = f(y_n + \Delta tk_3)\\
\end{align}
$$
