---
layout: single
title: "Taylor Series Expansion of Pressure"
permalink: /notes/FWguide/pressure/
author_profile: false
toc: true  # <--- Add this
---

<div class="toc-sidebar">
  <h3><i class="fa fa-bars"></i> Table of Contents</h3>
  <ul>
    <li><a href="#introduction-and-motivation">Introduction and Motivation</a></li>
  </ul>
</div>

<!-- NAVIGATION -->
<div class="nav-box-container">
  <div class="back-box">
    <div class="back-label"><i class="fa fa-arrow-left"></i> Back to:</div>
    <div class="back-link"><a href="/notes/FWguide/taylor/">Taylor Series Expansion of Velocities</a></div>
  </div>

  <div class="center-box">
    <a href="/notes/FWguide/intro/">FUNWAVE-TVD Main Page</a>
  </div>

  <div class="forward-box">
    <div class="forward-label">Forward to: <i class="fa fa-arrow-right"></i></div>
    <div class="forward-link"><a href="/notes/FWguide/intro/">Taylor Series Expansion of Pressure</a></div>
  </div>
</div>


## Introduction and Motivation

Conceptually, the goal here is to find an expression for pressure $p$ up to $\mathcal{O}(\mu^2)$ from the Taylor Series expansions discussed thus far, much like the last section. However, it turns out that the equations do not simplify very neatly, and require a fair amount of strategy in rearranging terms and making substitutions. 

It's first worth thinking about how pressure in such a model compares to linear (Airy) theory. In linear theory, the pressure is generally assumed to be hydrostatic such that $p=\rho g (\eta-z)$. This is exactly the weight of water above a given point in the water column. Pressure isn't even in the governing equations needed to solve for Airy Theory if we disregard pressure forcing, since the differential equation solved is just the Laplace equation for potential ($\nabla \phi = 0$) and none of the boundary conditions. 

Another issue with pressure is the fact that it's not associated with any time derivatives. In the Navier Stokes Equations, we only see it as spatial gradients. This means that it's much harder to apply algorithms that "step forward" in time (such as Euler's Method, Runge-Kutta, etc.). Often we have to approach solving the NSE in two stages, first where we step forward the velocities in time and then solve quite a challenging PDE to find the pressure at the time step. NSE solvers can spend upwards of 75% of their computational time on the pressure alone! Ideally, we are aiming to avoid such a computationally expensive option. Boussinesq theory allows us to do that with its Taylor Series truncated variables by finding an expression for pressure in terms of other known variables up to a given order.

## Depth Integration Setup
As is typical in our approach thus-far, we leverage the power of depth integration to get our pressure expression. The most natural thing to depth integrate is the $w$ term of the momentum balance, since a $\frac{\partial\mathbf{u}}{\partial z}$ appears there, and the integral of a derivative of $z$ has been our friend thus far since if $z$ is one of the bounds, we know that $p(z)$ will pop out. Also, the integral of $1$ is always quite easy too.

$$
\begin{aligned}
\underbrace{\epsilon \dfrac{\partial w}{\partial t} + \epsilon^2 u \dfrac{\partial w}{\partial x} + \epsilon^2 v \dfrac{\partial w}{\partial y} + \dfrac{\epsilon^2}{\mu^2} w \dfrac{\partial w}{\partial z}}_{\text{Not so friendly terms...}} + \underbrace{\epsilon \dfrac{\partial p}{\partial z} + 1 }_{\text{Friendly terms!}} &= 0
\end{aligned}
$$

However, the other remaining terms are decidedly *unfriendly*, being products of velocities and their derivatives from the Navier Stokes Equations. 

Another subtlety is that our choice in the bounds for the depth integration should probably be different. Recall how in the velocities, we integrated from $-h$ to $z$ and the kinematic bottom  boundary condition (KBBC) neatly popped up and reduced to 0. This was great for velocities, but pressure isn't even in the KBBC, so it's not as useful here. We *do* know pressure at the surface from the dynamic free surface boundary condition (DFSBC) though, which is at the surface! Thus, its more sensible to integrate from $z$ to $\epsilon\eta$, our scaled free surface. Doing this we have:

$$
\begin{aligned}
\underbrace{\epsilon \int_{z}^{\epsilon \eta}\dfrac{\partial w}{\partial t} \, dz'}_{(1)} + \underbrace{\epsilon^2 \int_{z}^{\epsilon \eta} ((\mathbf{u} \cdot \nabla_H)w)\, dz'}_{(2)} + \underbrace{\dfrac{\epsilon^2}{\mu^2} \int_{z}^{\epsilon \eta} \left(w \dfrac{\partial w}{\partial z'}\right)\, dz'}_{(3)} + \underbrace{\epsilon \int_{z}^{\epsilon \eta} \left(\dfrac{\partial p}{\partial z} \right)\, dz' + \int_{z}^{\epsilon \eta} 1 \, dz}_{(4)} &= 0
\end{aligned}
$$

We'll step through these in order of difficulty. Note that there a few different ways to rearrange terms to get to the same results.

## Integral 4 (Hydrostatic Pressure)

$$
\underbrace{\epsilon \int_{z}^{\epsilon \eta} \left(\dfrac{\partial p}{\partial z} \right)\, dz' + \int_{z}^{\epsilon \eta} 1 \, dz}_{(4)}
$$

This term contains the pressure gradient in the vertical as well as the effect of gravity, and is relatively straight forward to integrate:

$$
\begin{aligned}
\epsilon \int_{z}^{\epsilon \eta} \left(\dfrac{\partial p}{\partial z} \right)\, dz' + \int_{z}^{\epsilon \eta} 1 \, dz
&= \epsilon \int_{z}^{\epsilon \eta} \, dp + \int_{z}^{\epsilon \eta}\, dz' \\
&= \epsilon \int_{z}^{\epsilon \eta} \, dp + \int_{z}^{\epsilon \eta}\, dz'\\
&= \epsilon \left(\cancelto{0}{p(\epsilon \eta)} - p(z)\right) + \epsilon \eta - z
\end{aligned}
$$

This is where the DFBSC can be applied, where $p(\epsilon\eta)=0$. Thus we just have:

<div class="equation-block">
  <strong>Integral 4:</strong><br>
  $$
  (\epsilon \eta - z) - \epsilon p
  $$
</div>

The first part, $(\epsilon \eta - z)=p_h$ is exactly the hydrostatic pressure, consistent with Airy theory. This is just the weight of water above a point in the water column. The second term $\epsilon p$, is the (scaled) total pressure $p_{tot}$. In Airy Theory, $p_{tot} = p_h$. But here, we have those 3 other integrals that contribute terms. The combination of integrals 1,2,3 and thus the *hydrodynamic pressure*:

$$
p_{tot} = p_{h} + p_{d}
$$

Thus, solving through integrals 1-3 is effectively finding the hydrodynamic pressure.
## Integral 1 (Time Derivatives)

$$
\underbrace{\epsilon \int_{z}^{\epsilon \eta}\dfrac{\partial w}{\partial t} \, dz'}_{(1)} 
$$

There's a few ways we can go about rearranging terms here. While we *could* use the Liebnitz rule to take the derivative out of the integral, this is a bit uneccessary since we already have a nice solution for $w$ from the previous section. We can just differentiate that directly and plug into the integral. 


<div class="warning-block">
<strong>⚠ From here on out, time derivatives are denoted by a subscript t</strong><br>
</div>
Time derivatives can be brought inside spatial derivatives, and $z$ is unaffected by $\nabla_H$ or time derivatives and can be moved anywhere. Thus we have:
$$
\begin{aligned}
\frac{\partial w}{\partial t} = &w_t = \frac{\partial}{\partial t} \left[ -\mu^2 \nabla_H \cdot \left[\mathbf{u}_\alpha(z+h)\right]\right] \\
&w_t = \left[ -\mu^2 \nabla_H \cdot \left[\mathbf{u}_{\alpha t}(z+h)\right]\right] \\
&w_t = -\mu^2 z' (\nabla_H \cdot \mathbf{u}_{\alpha t}) - \mu^2[\nabla_H \cdot (h\mathbf{u}_{\alpha t})]\\
\end{aligned}
$$

Now we can integrate. As has been common, there is actually very little $z'$ dependency in the integrand, so large parts can be moved out

$$
\begin{aligned}
\epsilon \int_{z}^{\epsilon \eta}\dfrac{\partial w}{\partial t} \, dz' &= \epsilon \int_{z}^{\epsilon \eta} \left( -\mu^2 z (\nabla_H \cdot \mathbf{u}_{\alpha t}) - \mu^2[\nabla_H \cdot (h\mathbf{u}_{\alpha t})]\right) \\
&= -\mu^2\epsilon (\nabla_H \cdot \mathbf{u}_{\alpha t}) \int_{z}^{\epsilon \eta} z' \, dz'  -\mu^2\epsilon[\nabla_H \cdot (h\mathbf{u}_{\alpha t})]\int_{z}^{\epsilon \eta}\, dz'\\
\end{aligned}
$$

Thus we have:
<div class="equation-block">
  <strong>Integral 1:</strong><br>
  $$
  -\mu^2\epsilon (\nabla_H \cdot \mathbf{u}_{\alpha t}) \left(\frac{(\epsilon\eta)^2}{2} - \frac{z^2}{2}\right) -\mu^2\epsilon[\nabla_H \cdot (h\mathbf{u}_{\alpha t})](\epsilon\eta -z)\\
  $$
</div>

This integral contains all time derivatives within the hydrodynamic pressure.

## Integral 2 (Horizontal Advection)

$$
\underbrace{\epsilon^2 \int_{z}^{\epsilon \eta} ((\mathbf{u} \cdot \nabla_H)w)\, dz'}_{(2)}
$$

This is arguably the most complex integral, so we're going to do our best to analyze as little of it as possible explicitly. We do know $\mathbf{u}$ and $w$ and could plug in, but that gets incredibly messy. Instead, let's first note the order of the terms:

$$
\begin{aligned}
\mathbf{u} &= \mathbf{u}_\alpha + \mathcal{O}(\mu^2) \\
w &= \mathcal{O}(\mu^2) \\
\end{aligned}
$$

Since this term involves a product of $\mathbf{u}$ and $w$, the product of the $\mathcal{O}(\mu^2)$ part of $\mathbf{u}$ and $w$ become higher order terms, and are thus dropped. We see this as:


$$
(\mathbf{u} \cdot \nabla_H)w =  (\mathbf{u}_\alpha  \cdot \nabla_H)w + \mathcal{O}(\mu^4)
$$

We substitute this into the integrand now. Note that it's is a bit easier to work through the result term by term rather than in the $\nabla_H$ advective form, expanding  $\mathbf{u}$ to $(u_\alpha,v_\alpha)$ and considering each derivative separately:

$$
\int_{z}^{\epsilon \eta}(\mathbf{u}_\alpha  \cdot \nabla_H)w \, dz' = \int_{z}^{\epsilon \eta} \left( u_\alpha\frac{\partial w}{\partial x} + v_\alpha\frac{\partial w}{\partial y} \right)\, dz'
$$

We can substitute in $w$ to the derivative terms:

$$
\begin{aligned}
\frac{\partial w}{\partial x} &= -\mu^2\left[ z\frac{\partial}{\partial x}\left(\nabla_H \cdot u_\alpha\right) + \frac{\partial}{\partial x} \left(\nabla_H \cdot \left[hu_\alpha\right]\right)\right] \\
\frac{\partial w}{\partial y} &= -\mu^2\left[ z\frac{\partial}{\partial y}\left(\nabla_H \cdot v_\alpha\right) + \frac{\partial}{\partial y} \left(\nabla_H \cdot \left[hv_\alpha\right]\right)\right] \\
\end{aligned}
$$

Thus,

$$
\begin{aligned}
(\mathbf{u}_\alpha  \cdot \nabla_H)w &= u_\alpha\left\{ -\mu^2\left[ z\frac{\partial}{\partial x}\left(\nabla_H \cdot u_\alpha\right) + \frac{\partial}{\partial x} \left(\nabla_H \cdot \left[hu_\alpha\right]\right)\right]\right\} + \\
&= v_\alpha\left\{ -\mu^2\left[ z\frac{\partial}{\partial y}\left(\nabla_H \cdot v_\alpha\right) + \frac{\partial}{\partial y} \left(\nabla_H \cdot \left[hv_\alpha\right]\right)\right]\right\} \\
\end{aligned}
$$

Now we can depth-integrate.

$$
\begin{aligned}
\int_{z}^{\epsilon \eta} (\mathbf{u} \cdot \nabla_H) w\, dz' 
&= -\mu^2 u_\alpha \left\{ \frac{\partial}{\partial x} \left( \nabla \cdot u_\alpha \right) \left[ \int_{z}^{\epsilon \eta} z' dz' \right] + \frac{\partial}{\partial x} \left( \nabla \cdot (h u_\alpha) \right) \left[ \int_{z}^{\epsilon \eta} dz' \right] \right\} \\
&\quad - \mu^2 v_\alpha \left\{ \frac{\partial}{\partial y} \left( \nabla \cdot v_\alpha \right) \left[ \int_{z}^{\epsilon \eta} z' dz' \right] + \frac{\partial}{\partial y} \left( \nabla \cdot (h v_\alpha) \right) \left[ \int_{z}^{\epsilon \eta} dz' \right] \right\} \\
&= -\mu^2 u_\alpha \left\{ \frac{\partial}{\partial x} \left( \nabla \cdot u_\alpha \right) \left[ \frac{(\epsilon \eta)^2 - z^2}{2} \right] + \frac{\partial}{\partial x} \left( \nabla \cdot (h u_\alpha) \right) [\epsilon \eta - z] \right\} \\
&\quad -\mu^2 v_\alpha \left\{ \frac{\partial}{\partial y} \left( \nabla \cdot v_\alpha \right) \left[ \frac{(\epsilon \eta)^2 - z^2}{2} \right] + \frac{\partial}{\partial y} \left( \nabla \cdot (h v_\alpha) \right) [\epsilon \eta - z] \right\}
\end{aligned}
$$

Condensing this as much as possible and evaluating the integrals, we have:


<div class="equation-block">
  <strong>Integral 2:</strong><br>
$$
-\mu^2\epsilon^2\mathbf{u}_\alpha \cdot \nabla_H \left( \nabla_H \cdot \mathbf{u}_\alpha \right) \left( \frac{(\epsilon\eta)^2 - z^2}{2} \right) 
-\mu^2\epsilon^2 \mathbf{u}_\alpha \cdot \nabla_H \left( \nabla_H \cdot (h \mathbf{u}_\alpha) \right) (\epsilon\eta - z)
$$
</div>

## Integral 3 (Vertical Advection)

$$
\underbrace{\dfrac{\epsilon^2}{\mu^2} \int_{z}^{\epsilon \eta} \left(w \dfrac{\partial w}{\partial z'}\right)\, dz'}_{(3)}
$$

This is a bit simpler. We'll start by finding $\frac{\partial w}{\partial z}$ and then $w\frac{\partial w}{\partial z}$

$$
\begin{aligned}
\frac{\partial w}{\partial z} &= \frac{\partial}{\partial z}\left[ -\mu^2 \nabla_H \cdot \left[\mathbf{u}_\alpha(z+h)\right]\right] = -\mu^2 \left[\nabla_H \cdot \mathbf{u}_\alpha\right] \\
w\frac{\partial w}{\partial z} &= \left[ -\mu^2 \nabla_H \cdot \left[\mathbf{u}_\alpha(z+h)\right]\right] \left(-\mu^2 \left[\nabla_H \cdot \mathbf{u}_\alpha\right]\right) \\
&= \mu^4 \left\{ z \left[ \nabla_H \cdot \mathbf{u}_\alpha \right]^2 
+ \left[ \nabla_H \cdot (h \mathbf{u}_\alpha) \right] \left[ \nabla_H \cdot (\mathbf{u}_\alpha) \right] \right\}
\end{aligned}
$$

The fact that this is $\mathcal{O}(\mu^4)$ is okay here since in the integral we multiply by $\frac{\epsilon^2}{\mu^2}$. Thus, we continue on with depth integration:

$$
\begin{aligned}
\int_{z}^{\epsilon \eta} \left[ w \frac{\partial w}{\partial z'} \right] dz'
&= \mu^4 \left\{ \left[ \nabla_H \cdot \mathbf{u}_\alpha \right]^2 \int_{z}^{\epsilon \eta} z' \, dz' 
+ \left[ \nabla_H \cdot (h \mathbf{u}_\alpha) \right] \left[ \nabla_H \cdot (\mathbf{u}_\alpha) \right] \int_{z}^{\epsilon \eta} dz' \right\} \\
&= \mu^4 \left\{ \left[ \nabla_H \cdot \mathbf{u}_\alpha \right]^2 \left( \frac{\eta^2 - z^2}{2} \right)
+ \left[ \nabla_H \cdot (h \mathbf{u}_\alpha) \right] \left[ \nabla_H \cdot (\mathbf{u}_\alpha) \right] (\eta - z) \right\}
\end{aligned}
$$

<div class="equation-block">
  <strong>Integral 3:</strong><br>
$$
\dfrac{\epsilon^2}{\mu^2} \int_{z}^{\epsilon \eta} \left(w \dfrac{\partial w}{\partial z'}\right)\, dz' = \epsilon^2\mu^2 \left[ \nabla_H \cdot \mathbf{u}_\alpha \right]^2 \left( \frac{\eta^2 - z^2}{2} \right)
+ \epsilon^2\mu^2\left[ \nabla_H \cdot (h \mathbf{u}_\alpha) \right] \left[ \nabla_H \cdot (\mathbf{u}_\alpha) \right] (\eta - z) 
$$
</div>

## Assembling the Terms
Now, with all the integrals calculated, we can assemble them to find the total pressure. From Term (4) we have pressure $p$ on the LHS with the hydrostatic pressure:

$$
\begin{aligned}
(\epsilon \eta - z) - \epsilon p &=  -\mu^2\epsilon (\nabla_H \cdot \mathbf{u}_{\alpha t}) \left(\frac{(\epsilon\eta)^2}{2} \\- \frac{z^2}{2}\right) -\mu^2\epsilon[\nabla_H \cdot (h\mathbf{u}_{\alpha t})](\epsilon\eta -z)\\
  &  -\mu^2\epsilon^2\mathbf{u}_\alpha \cdot \nabla_H \left( \nabla_H \cdot \mathbf{u}_\alpha \right) \left( \frac{(\epsilon\eta)^2 - z^2}{2} \right) 
-\mu^2\epsilon^2 \mathbf{u}_\alpha \cdot \nabla_H \left( \nabla_H \cdot (h \mathbf{u}_\alpha) \right) (\epsilon\eta - z) \\
& +\epsilon^2\mu^2 \left[ \nabla_H \cdot \mathbf{u}_\alpha \right]^2 \left( \frac{\eta^2 - z^2}{2} \right)
+ \left[ \nabla_H \cdot (h \mathbf{u}_\alpha) \right] \left[ \nabla_H \cdot (\mathbf{u}_\alpha) \right] (\eta - z) 
\end{aligned}
$$

Now we just rearrange, very carefully. We start by getting pressure alone and dividing through by $\epsilon$ and cleaning up the $\frac{1}{2}$ terms and factoring out a negative from them

$$
\begin{aligned}
p &= \left(\eta - \frac{z}{\epsilon}\right) -\mu^2(\nabla_H \cdot \mathbf{u}_{\alpha t}) \left(\tfrac{1}{2}\left[z^2 -(\epsilon\eta)^2  \right]\right) -\mu^2[\nabla_H \cdot (h\mathbf{u}_{\alpha t})](z-\epsilon \eta)\\
  &  -\mu^2\epsilon\mathbf{u}_\alpha \cdot \nabla_H \left( \nabla_H \cdot \mathbf{u}_\alpha \right) \left(\tfrac{1}{2}\left[z^2 -(\epsilon\eta)^2  \right]\right) 
-\mu^2\epsilon \mathbf{u}_\alpha \cdot \nabla_H \left( \nabla_H \cdot (h \mathbf{u}_\alpha) \right) (z-\epsilon \eta) \\
& +\epsilon\mu^2 \left[ \nabla_H \cdot \mathbf{u}_\alpha \right]^2 \left(\tfrac{1}{2}\left[z^2 -(\epsilon\eta)^2  \right]\right)
- \epsilon\mu^2\left[ \nabla_H \cdot (h \mathbf{u}_\alpha) \right] \left[ \nabla_H \cdot (\mathbf{u}_\alpha) \right] (z-\epsilon \eta) 
\end{aligned}
$$

Next, factor out the $\mu^2$:

$$
\begin{aligned}
p = \left(\eta - \frac{z}{\epsilon}\right) -&\mu^2\left\{  (\nabla_H \cdot \mathbf{u}_{\alpha t}) \left(\tfrac{1}{2}\left[z^2 -(\epsilon\eta)^2  \right]\right) \right.\\
&+ [\nabla_H \cdot (h\mathbf{u}_{\alpha t})](z-\epsilon \eta)\\
  &  + \epsilon\mathbf{u}_\alpha \cdot \nabla_H \left( \nabla_H \cdot \mathbf{u}_\alpha \right) \left(\tfrac{1}{2}\left[z^2 -(\epsilon\eta)^2  \right]\right) \\
& +\epsilon \mathbf{u}_\alpha \cdot \nabla_H \left( \nabla_H \cdot (h \mathbf{u}_\alpha) \right) (z-\epsilon \eta) \\
& -\epsilon \left[ \nabla_H \cdot \mathbf{u}_\alpha \right]^2 \left(\tfrac{1}{2}\left[z^2 -(\epsilon\eta)^2  \right]\right)\\
& + \left.\epsilon\left[ \nabla_H \cdot (h \mathbf{u}_\alpha) \right] \left[ \nabla_H \cdot (\mathbf{u}_\alpha) \right] (z-\epsilon \eta) \right\}
\end{aligned}
$$

Shuffling around some terms even more, we come to 



<!-- NAVIGATION -->
<div class="nav-box-container">
  <div class="back-box">
    <div class="back-label"><i class="fa fa-arrow-left"></i> Back to:</div>
    <div class="back-link"><a href="/notes/FWguide/intro/">Taylor Series Expansion of Pressure</a></div>
  </div>

  <div class="center-box">
    <a href="/notes/FWguide/intro/">FUNWAVE-TVD Main Page</a>
  </div>

  <div class="forward-box">
    <div class="forward-label">Forward to: <i class="fa fa-arrow-right"></i></div>
    <div class="forward-link"><a href="/notes/FWguide/next/">Taylor Series Expansion of Pressure</a></div>
  </div>
</div>
