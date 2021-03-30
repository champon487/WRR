# 1D Bed evolution model

## Flow model

The shallow water flow model (St. Venant) is used to calculate the 1D flow field. The governing equation is descritized on Staggered grid system.

### Governing equation

$$
\frac{\partial h}{\partial t}+\frac{\partial uh}{\partial x} = 0
$$

![\frac{\partial h}{\partial t}+\frac{\partial uh}{\partial x} = 0
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cfrac%7B%5Cpartial+h%7D%7B%5Cpartial+t%7D%2B%5Cfrac%7B%5Cpartial+uh%7D%7B%5Cpartial+x%7D+%3D+0%0A)

$$
\frac{\partial u}{\partial t}+u\frac{\partial u}{\partial x} = -g\frac{\partial H}{\partial x} - \frac{gn^2u|u|}{h^{1/3}}
$$
![\frac{\partial u}{\partial t}+u\frac{\partial u}{\partial x} = -g\frac{\partial H}{\partial x} - \frac{gn^2u|u|}{h^{1/3}}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cfrac%7B%5Cpartial+u%7D%7B%5Cpartial+t%7D%2Bu%5Cfrac%7B%5Cpartial+u%7D%7B%5Cpartial+x%7D+%3D+-g%5Cfrac%7B%5Cpartial+H%7D%7B%5Cpartial+x%7D+-+%5Cfrac%7Bgn%5E2u%7Cu%7C%7D%7Bh%5E%7B1%2F3%7D%7D%0A)

where, t: time, x: downstream cordinate, u: flow velocity, h: water depth, H: water surface elevation, g: gravitational acceleration, and n: Manning's roughness coefficient.

## Morphodynamics

The model deals with bed- and suspended-load transport. The bed evolutaiton is then calculated using Exner equation.

### Governing equation

$$
(1-\lambda)\frac{\partial \eta}{\partial t}+\frac{\partial q_b}{\partial x} = D-E
$$

$$
q_b = 4(\theta-\theta_c)^{3/2}\sqrt{Rgd^3}
$$

where, $\eta$: bed elevation, $\theta$: Shields number, $\theta_c$: criticla Shields number, $R$: specific weight of sediment in fluid, $d$: grain size, $q_b$: bedload flux, $E$: entrainment rate, $D$: deposition rate.