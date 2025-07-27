---
layout: single
title: "The Momentum Balance Equation (w/o Vorticity)"
permalink: /notes/FWguide/momentum/
author_profile: false
toc: true  # <--- Add this
---


## Simplifying Notation

Due to the extensive length of the terms, it is useful to introduce some simplifying notation to group like terms. We define some notation for common combinations of variables to simplify expressions. In introducing such notation, it is important to remember which terms contain spatial and time
dependencies.



<div class="equation-block">
  <strong>Definition of $A$ and $B$</strong><br>
$$
\begin{aligned}
A &\equiv \nabla_H \cdot (h \mathbf{u_\alpha}) \\
B &\equiv \nabla_H \cdot \mathbf{u_\alpha}
\end{aligned}
$$

where both $A$ and $B$ contain a time dependency due to $\mathbf{u_\alpha}$, and a spatial dependency due to both $h$ and $\mathbf{u_\alpha}$
</div>

In terms of $A$ and $B$ then, the hydrodynamic variables are:

<div class="equation-block">
  <strong>Velocities in terms of $A$ and $B$</strong><br>
$$
\begin{aligned}
\mathbf{u} &= \mathbf{u_\alpha} + \mu^2(z_\alpha-z)\nabla_H A + \mu^2\tfrac{1}{2}(z_\alpha^2-z^2)\nabla_H B\\
w &= -\mu^2(A+zB) \\
\end{aligned}
$$
<strong>Pressure in terms of $A$ and $B$</strong><br>
$$
\begin{aligned}
p = g(\eta-z) + (z-\eta) & \left\{A_t + \tfrac{1}{2}(z+\eta)B_t   \right. \\
& \quad \mathbf{u}_\alpha \cdot \nabla_H A + \tfrac{1}{2}(z+\eta)\mathbf{u}_\alpha \cdot \nabla_H B - \\
& \quad \left. \left[ AB - \tfrac{1}{2}(z+\eta)B^2\right]\right\}
\end{aligned} 
$$

where both $A$ and $B$ contain a time dependency due to $\mathbf{u_\alpha}$, and a spatial dependency due to both $h$ and $\mathbf{u_\alpha}$
</div>


## Substitution into Euler Equations
With all the variables solved for and simplified, all that remains is plugging them into the Euler Equations. Here, we will use just the horizontal Euler Equations and solve term by term as indicated. Note that the vorticity term $\mathcal{V}$ is quite complicated and will be addressed later.


$$
\underbrace{\frac{\partial \mathbf{u}}{\partial t}}_{(1)} + \underbrace{\nabla_H (\mathbf{u} \cdot \mathbf{u})}_{(2)} + \underbrace{\frac{\partial}{\partial z} (\tfrac{1}{2}w \cdot w)}_{(3)} + \underbrace{\mathbf{\Omega} \times \mathbf{q}}_{\mathcal{V}} + \underbrace{\tfrac{1}{2}\nabla_H p}_{(4)} = 0 
$$

The general strategy within each term will be to separate out/group the $z_\alpha$ and the $z$ terms.
Additional notation will be introduced for each group formed.

## Term 1

