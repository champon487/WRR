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


## Example1: 1D bed evolution after a cutoff

Discharge       = 200 m<sup>3</sup>  
Channel width   = 200 m  
Manning's n     = 0.03  
Bed slope       = 0.002 (mild reach), 0.01 (steep reach)  
Grain size      = 2 mm  

[![](https://img.youtube.com/vi/W-FToCP7g84/0.jpg)](https://www.youtube.com/watch?v=W-FToCP7g84)

## Example2: Sedimentation in a dam

Discharge       = 200 m<sup>3</sup>  
Channel width   = 200 m  
Manning's n     = 0.03  
Bed slope       = 0.002  
Grain size      = 0.2 mm  
Downstream WL   = 3 m (fixed)

[![](https://img.youtube.com/vi/q0Y_GjS6tMg/0.jpg)](https://www.youtube.com/watch?v=q0Y_GjS6tMg)