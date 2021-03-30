# 1D Bed evolution model

## Flow model

The shallow water flow model (St. Venant) is used to calculate the 1D flow field. The governing equation is descritized on Staggered grid system.

### Governing equation

<img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;h}{\partial&space;t}&plus;\frac{\partial&space;uh}{\partial&space;x}&space;=&space;0" />

<img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;u}{\partial&space;t}&plus;u\frac{\partial&space;u}{\partial&space;x}&space;=&space;-g\frac{\partial&space;H}{\partial&space;x}&space;-&space;\frac{gn^2u|u|}{h^{1/3}}" />

where, t: time, x: downstream cordinate, u: flow velocity, h: water depth, H: water surface elevation, g: gravitational acceleration, and n: Manning's roughness coefficient.

## Morphodynamics

The model deals with bed- and suspended-load transport. The bed evolutaiton is then calculated using Exner equation.

### Governing equation

<img src="https://latex.codecogs.com/gif.latex?(1-\lambda)\frac{\partial&space;\eta}{\partial&space;t}&plus;\frac{\partial&space;q_b}{\partial&space;x}&space;=&space;D-E" />

<img src="https://latex.codecogs.com/gif.latex?q_b&space;=&space;4(\theta-\theta_c)^{3/2}\sqrt{Rgd^3}" />

where, &eta;: bed elevation, &theta;: Shields number, &theta;<sub>c</sub>;: criticla Shields number, R: specific weight of sediment in fluid, d: grain size, q<sub>b</sub>: bedload flux, E: entrainment rate, D: deposition rate.
