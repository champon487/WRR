# 1D Bed evolution model

## Flow model

The shallow water flow model (St. Venant) is used to calculate the 1D flow field. The governing equation is descritized on Staggered grid system.

### Governing equation

$$
\frac{\partial h}{\partial t}+\frac{\partial uh}{\partial x} = 0
$$

<img src="https://latex.codecogs.com/gif.latex?\int_a^bf(x)dx" />

<img src="https://latex.codecogs.com/gif.latex?\frac{\partialh}{\partialt}+\frac{\partialuh}{\partialx}=0" />

$$
\frac{\partial u}{\partial t}+u\frac{\partial u}{\partial x} = -g\frac{\partial H}{\partial x} - \frac{gn^2u|u|}{h^{1/3}}
$$
<img src="https://latex.codecogs.com/gif.latex?\frac{\partial u}{\partial t}+u\frac{\partial u}{\partial x} = -g\frac{\partial H}{\partial x} - \frac{gn^2u|u|}{h^{1/3}}" />

where, t: time, x: downstream cordinate, u: flow velocity, h: water depth, H: water surface elevation, g: gravitational acceleration, and n: Manning's roughness coefficient.

## Morphodynamics

The model deals with bed- and suspended-load transport. The bed evolutaiton is then calculated using Exner equation.

### Governing equation

$$
(1-\lambda)\frac{\partial \eta}{\partial t}+\frac{\partial q_b}{\partial x} = D-E
$$

<img src="https://latex.codecogs.com/gif.latex?(1-\lambda)\frac{\partial \eta}{\partial t}+\frac{\partial q_b}{\partial x} = D-E" />

$$
q_b = 4(\theta-\theta_c)^{3/2}\sqrt{Rgd^3}
$$

<img src="https://latex.codecogs.com/gif.latex?q_b = 4(\theta-\theta_c)^{3/2}\sqrt{Rgd^3}" />

where, $\eta$: bed elevation, $\theta$: Shields number, $\theta_c$: criticla Shields number, $R$: specific weight of sediment in fluid, $d$: grain size, $q_b$: bedload flux, $E$: entrainment rate, $D$: deposition rate.
