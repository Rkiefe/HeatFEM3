# HeatFEM
HeatFEM is an open-source implementation of the finite element method for 3D objects. Users can import stl files to simulate complex shapes, set the initial temperature profile, thermal conductivity, time of the simulation and more, to simulate heat on complex 3D shapes.

Validation of HeatFEM was performed by comparing with an analytic solution,
$T = T_0 e^{-t} cos(x+y+z)$
with boundary conditions:
$-k\nabla T \cdot \vec{n} = kT_0 e^{-t}sin(x+y+z) [1,1,1]\cdot \vec{n}$

<img src="https://github.com/user-attachments/assets/824870fa-2b5a-4a0b-a339-ce0ca52be7b5" alt="FEM" width="400" height="300"/>
<img src="https://github.com/user-attachments/assets/ea17495c-0b04-423d-8281-8de3214d2ae9" alt="FEM" width="400" height="300"/>

## Model
The heat equation is as follows,

$\rho C_p \frac{\partial T}{\partial t} - \nabla \cdot (k\nabla T) = \frac{d q_V}{dt}$

With Neumann boundary conditions: 

$k\nabla T \cdot \vec{n} = h\left(T_{ext}-T\right)$
