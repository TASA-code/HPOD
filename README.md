# Resulting Plot

<p align="center">
  <img src="plot/3d.png" width="900" title="C++">
</p>


# Orbit_Determination


This code serves the purpose for orbit determination of a satellite orbiting aroung Earth with the effect of J2 perturbation. 
The code takes in classical orbital element (SMA, e, i, RAAN, &omega;, &nu;) into account to calculate and iterate the final position and velocity vector

- SMA  : Semi-Major axis,
- e    : Eccentricity,
- i    : Inclination,
- RAAN : Right Ascension of Ascending Node,
- &omega;: Argument of periapsi,
- &nu; : True Anomaly.

This code takes in the argument from the file "input.txt" in the root folder
> [!IMPORTANT]
> Specific format must be following as indicated in the example input.txt file

# Usage
The repository compiles using cMake. To compile the program, move to the build directory

```command line
mkdir build
cd build
```

and run 
```command line
cmake ..
make test
```
The target test case will propagte the orbit as well as plot and open the plot for visualisation when completed.



# Mathematical Formulas
## Orbital Equation

In this section, we provide two important formulas related to orbital mechanics.

The formula for calculating the orbital radius (r) is given by:

$$
r = \frac{h^2}{\mu}\frac{1}{1+e\cos\theta}(\cos\theta\hspace{0.2cm}\mathbf{i}_e + \sin\theta\hspace{0.2cm}\mathbf{i}_p)
$$

The formula for calculating the orbital velcotiy (v) is given by:

$$
v = \frac{\mu}{h}(\sin\theta\hspace{0.2cm}\mathbf{i}_e + (e+\cos\theta)\hspace{0.2cm}\mathbf{i}_p)
$$

> [!IMPORTANT]
> The orbital coordinate are transformed from Perifocal to ECI

## Gravitational Acceleration and Oblateness Perturbation (J2)

The orbital equation under gravitational acceleration and oblateness perturbation shown below,

$$
\frac{d^2r}{dt^2} + \mu\frac{\mathbf{r}}{r^3} = \mathbf{a}_d = -\frac{3}{2}\frac{J_2\mu R_E\mathbf{r}}{2r^5} = \begin{bmatrix} 1-\frac{5r_z^2}{r^2} \\\ 1-\frac{5r_z^2}{r^2} \\\ 3-\frac{5r_z^2}{r^2} \end{bmatrix}
$$


# Time Integration

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


