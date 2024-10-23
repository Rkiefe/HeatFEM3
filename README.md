# HeatFEM
This is a finite element implementation of the 3D heat equation, with Neumann boundary conditions.
![example](https://github.com/user-attachments/assets/66892836-38cd-4071-a49e-d64b8f8df0c0)

## Model
The heat equation is as follows,

$\rho C_p \frac{\partial T}{\partial t} - \nabla \cdot (k\nabla T) = \frac{d q_V}{dt}$

With Neumann boundary conditions: 

$k\nabla T \cdot \vec{n} = h\left(T_{ext}-T\right)$
