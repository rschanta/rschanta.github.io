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
  <strong>Definition of A and B</strong><br>
$$
\begin{aligned}
A &\equiv \nabla_H \cdot (h \mathbf{u_\alpha}) \\
B &\equiv \nabla_H \cdot \mathbf{u_\alpha}
\end{aligned}
$$

where both $A$ and $B$ contain a time dependency due to $\mathbf{u_\alpha}$, and a spatial dependency due to both $h$ and $\mathbf{u_\alpha}$
</div>

## Substitution into Euler Equations