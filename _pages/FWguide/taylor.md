---
layout: single
title: "Taylor Series Expansion of Velocities"
permalink: /notes/FWguide/taylor/
author_profile: false
toc: true  # <--- Add this
---
<div class="toc-sidebar">
  <h3><i class="fa fa-bars"></i> Table of Contents</h3>
  <ul>
    <li><a href="#motivation">Motivation</a></li>
    <li><a href="#goals-and-strategy">Goals and Strategy</a></li>
    <li><a href="#location-of-z_alpha-and-notation">Location of $z_\alpha$ and Notation</a></li>
    <li><a href="#horizontal-vorticity-assumption">Horizontal Vorticity Assumption</a></li>
    <li><a href="#relating-w-to-u">Relating $w(z)$ to $\mathbf{u}$</a></li>
    <li><a href="#taylor-series-expansions-of-velocities">Taylor Series Expansions of Velocities</a></li>
    <li><a href="#extension-to-an-arbitrary-depth">Extension to an Arbitrary Depth</a></li>
    <li><a href="#implications-of-the-taylor-series-expansion">Implications</a></li>
    <li><a href="#conclusion">Conclusion</a></li>
  </ul>
</div>


<!-- NAVIGATION -->
<div class="nav-box-container">
  <div class="back-box">
    <div class="back-label"><i class="fa fa-arrow-left"></i> Back to:</div>
    <div class="back-link"><a href="/notes/FWguide/intro/">Overview Page</a></div>
  </div>

  <div class="center-box">
    <a href="/notes/FWguide/intro/">FUNWAVE-TVD Main Page</a>
  </div>

  <div class="forward-box">
    <div class="forward-label">Forward to: <i class="fa fa-arrow-right"></i></div>
    <div class="forward-link"><a href="/notes/FWguide/pressure/">Taylor Series Expansion of Pressure</a></div>
  </div>
</div>


## Motivation
So far, all of our equations are in three spatial dimensions- $(x,y,z)$ as well as time $t$. When discussing the fluid velocities $(u,v,w)$ and pressure $p$, this is sensible, as we expect these hydrodynamic variables to vary in the horizontal as waves propagate and vary in the vertical in response to pressure, boundary conditions, etc. 

Yet the variable $\eta$ is inherently different, as it is really a surface in $z$. We have the surface $z=\eta(x,y,t)$ representing the top of the water column. We inherently think of wave propagation in a 2D plane. Following linear theory, we also assumed a scale of $k_0$ for $(x,y)$ and $h_0$ for $z$. So there is a very natural split between the horizontal and vertical dimensions of the problem.

This split is one such motivation in the formulation Boussinesq models. If our goal is to model 2D wave propagation, why not simplify the $z$ dimension? Historically, *depth integration* has been a popular tool for simplifying the $z$ dimension (see), but Boussinesq models instead use a **Taylor Series** in the $z$ dimension about some reference elevation in the water column denoted $z_\alpha$:



<div style="text-align: center;">
  <img src="{{ '/assets/images/FWguide/taylor.PNG' | relative_url }}" alt="Alt Text" width="60%">
</div>

In general, such a Taylor Series will have the form of:

<div class="equation-block">
  <strong>Taylor Series Expansion about $z_\alpha$:</strong><br>
  $$
  \begin{aligned}
  \mathbf{q} = \mathbf{q}(x,y,-z_\alpha,t) + (z_\alpha-z)\left(\frac{\partial \mathbf{q}}{\partial z}\right)_{z=-z_\alpha}+\frac{1}{2}(z_\alpha-z)^2\left(\frac{\partial^2 \mathbf{q}}{\partial z^2}\right)_{z=-z_\alpha} + ...\\

  \end{aligned}
  $$
</div>

It should be noted that this is only a Taylor Series expansion in $z$- we maintain the full dimensionality of $(x,y)$. As is typical in multivariable calculus, the other variables we're not working with are just ``along for the ride" so to speak. 

We can extend this Taylor Series as far as we'd like. The more terms we add, the more accurate our expansion will be, but it will also be more complicated and computationally expensive. In the case of FUNWAVE-TVD, we keep terms up to $\mathcal{O}(\mu^2)$. It will become clear that as we add higher and higher terms to a given variable, the $\mu$ parameter tends to increase with each order. If we were to keep terms only up to $\mathcal{O}(1)$, we would revert back to simply linear theory. Higher order Boussinesq models have been developed (See Gobbi), but these formulations tend to be extremely complicated and computationally very expensive so are infrequently used.


## Goals and Strategy

So the goal is to the take the governing equations below and use truncated Taylor Series of the variables $(u,v,w)$ and $p$. 

<div class="equation-block">
  <strong>Continuity Equation:</strong><br>
  $$
  \begin{aligned}
  \mu^2 \left( \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} + \frac{\partial w}{\partial z}\right) &= 0 \\

  \end{aligned}
  $$
  <strong>Momentum Equations:</strong><br>
  $$
  \begin{aligned}
  \mu^2 \dfrac{\partial u}{\partial t} + \epsilon \mu^2 u \dfrac{\partial u}{\partial x} + \epsilon \mu^2 v \dfrac{\partial u}{\partial y} + \epsilon w \dfrac{\partial u}{\partial z} + \mu^2 \dfrac{\partial p}{\partial x} &= 0 \\
  \mu^2 \dfrac{\partial v}{\partial t} + \epsilon \mu^2 u \dfrac{\partial v}{\partial x} + \epsilon \mu^2 v \dfrac{\partial v}{\partial y} + \epsilon w \dfrac{\partial v}{\partial z} + \mu^2 \dfrac{\partial p}{\partial y} &= 0 \\
  \epsilon \dfrac{\partial w}{\partial t} + \epsilon^2 u \dfrac{\partial w}{\partial x} + \epsilon^2 v \dfrac{\partial w}{\partial y} + \dfrac{\epsilon^2}{\mu^2} w \dfrac{\partial w}{\partial z} + \epsilon \dfrac{\partial p}{\partial z} + 1 &= 0
  \end{aligned}
  $$
</div>

 
This is easier said than done. One of the main challenges is the fact that we need expressions for the derivatives $\frac{\partial^n \mathbf{q}}{\partial z^n}$ at the reference elevation $z=-z_\alpha$, which we don't inherently have neat expressions for as-is. Also, our goal is a 2DH model, so we want to eliminate instances of the vertical $z$ and $w$ where they appear, and they appear in just about every equation right now. So it's largely a game of playing around with our 4 governing equations and boundary conditions to achieve these find these derivative terms, eliminate $z$ and $w$ terms, and cross out higher order terms along the way.

One tool that ends up being extremely useful is *depth-integration*, whereby we integrate some expression in the $z$ dimension:

<div class="equation-block">
  <strong>Depth Integration:</strong><br>
  $$
  \begin{aligned}
  \int_a^b f(x,y,z',t)dz'\\

  \end{aligned}
  $$
</div>

## Location of $z_\alpha$ and Notation
The specific choice of what exactly $z_\alpha$ should be varies throughout the literature. Historically, it has been common to expand the variables at $z_\alpha=-h$ such that the expansion is about the bed:

Expressions tend to get very long from here on out, so the following shorthands are adopted:

- The subscript $B$ is used for any variable evaluated at $z=-h$ but left as a function of $x,y,t$. For example, $u_B=u(x,y,z=-h,t)$. 
- The subscript $\alpha$ is used for any variable evaluated at $z=-z_\alpha$ but left as a function of $x,y,t$. For example, $u_\alpha=u(x,y,z=-z_\alpha,t)$. 
- Unless otherwise specified, any variable without arguments explicitly stated is assumed to be a function of $x,y,z,t$. For example,$u=u(x,y,z,t)$. 

⚠ In general, it is important to keep track of what exactly is a function of what and in what dimensions, so special care should be taken with shorthand notation!

## Horizontal Vorticity Assumption

We will begin by considering the Taylor Series expansion of the vertical velocity $\mathbf{u}=(u,v)$ about the bed. For now, we will consider just a first order expansion: 

$$
  \begin{aligned}
  \mathbf{u}= \mathbf{u}_B + (z+h)\left(\frac{\partial \mathbf{u}}{\partial\mathbf{z}}\right)_{B}+ ...\\

  \end{aligned}
  $$

$\mathbf{u}\_B$ is ultimately what we're trying to model, so the only unknown here is the derivative term $\left(\frac{\partial \mathbf{u}}{\partial\mathbf{z}}\right)_{B}$. This is a problem as-is, since we have no way of expressing this term. We need to invoke another assumption on vorticity given by:

<div class="equation-block">
  <strong>Zero Horizontal Vorticity Assumption</strong><br>
  The field is irrotational in the horizontal such that:
  $$
  \begin{aligned}
  \left[
      \left( \frac{\partial w}{\partial y} - \frac{\partial v}{\partial z} \right),
      \left( \frac{\partial u}{\partial z} - \frac{\partial w}{\partial x} \right),
      \left( \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y} \right)
      \right]
      =
      \left[
      0,\ 0,\ \left( \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y} \right)
      \right]\\

  \end{aligned}
  $$
</div>

This assumption has important consequences. Physically, this would imply that any rotational motion of the fluid occurs in the horizontal plane, spinning around a vertical axis. Vertical shear (such as $\frac{\partial u}{\partial z}$) is in general quite weak.

Importantly for our derivation, we now have a way to find an expression for $\left(\frac{\partial \mathbf{u}}{\partial\mathbf{z}}\right)_{B}$. Looking at the first two terms, it's clear that:

$$
  \begin{aligned}
      \frac{\partial v}{\partial z} = \frac{\partial w}{\partial y} \\
       \frac{\partial u}{\partial z} = \frac{\partial w}{\partial x}
  \end{aligned}
  $$

  Which in more compact vector notation is given as :

$$
\begin{aligned}
\frac{\partial \mathbf{u}}{\partial z} = \nabla_H w\\
\end{aligned}
$$

Now, we’ve related the vertical derivatives of the horizontal velocities to horizontal derivatives of the vertical velocity, which should be easily obtainable if we have an expression for $w(z)$. This is possible through depth integration.

## Relating $w(z)$ to $\mathbf{u}$
The most direct way to obtain some expression for $w(z)$ is through the depth-integrated contuinity equation. Let's condense the horizontal velocities into $\mathbf{u}$

$$
\begin{aligned}
\mu^2 \left( \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} \right) + \frac{\partial w}{\partial z} = \mu^2 (\nabla_H \cdot \mathbf{u}) + \frac{\partial w}{\partial z} = 0\\
\end{aligned}
$$

$w(z)$ appears as a derivative of z, so to isolate it, we need to depth-integrate the expression. At least one of the bounds of the depth integral needs to be $z$ itself if we want the output to also be a function of $z$. Thus, we use a dummy variable $z'$ within the integral itself. The second bound can essentially be whatever is convenient. For reasons that will become clear upon the integration, using $−h$ as a lower bound generally is a good choice:

$$
\begin{aligned}
\underbrace{\mu^2 \int_{-h}^{z} \mu^2 (\nabla_H \cdot \mathbf{u})\, dz'}_{(1)} + \underbrace{\int_{-h}^{z} \frac{\partial w}{\partial z'}\, dz'}_{(2)} = 0
\end{aligned}
$$

Term (1) is a classic use-case of Liebnitz Rule- we have the integral of a derivative (or a divergence here since $\textbf{u}$ is a vector). 

$$
\begin{aligned}
&\underbrace{\mu^2 \nabla_H \cdot \int_{-h}^{z} \mathbf{u} \, dz'}_{(1a)} 
+ \underbrace{\cancelto{0}{\mu^2 \mathbf{u}(z)(\nabla_H \cdot z)}}_{(1b)}  
+ \underbrace{\mu^2 \mathbf{u}(-h)(\nabla_H \cdot h) }_{(1c)} 
\end{aligned}
$$

Term (1b) is a bit strange- it is the derivative of an independent variable in space with respect
to another independent variable in space. If we assume that this $z$ value represents a specific plane (e.g. $z=-5$), then it doesn't vary in the $x-y$ plane and is thus $0$. So the term is removed.

Now if we turn to Term (2), this integral is much simpler since it reduces to a integral in $dw$:

$$
\begin{aligned}
\int_{-h}^{z} \frac{\partial w}{\partial z'}\, dz' = \int_{-h}^{z} dw = w(z) - w(-h)
\end{aligned}
$$

Now combining these simplified expressions, we have:

$$
\begin{aligned}
&\mu^2 \nabla_H \cdot \int_{-h}^{z} \mathbf{u} \, dz'
+ \mu^2 \mathbf{u}(-h)(\nabla_H \cdot h) + w(z) - w(-h) = 0
\end{aligned}
$$

If we rearrange just a bit, we notice that the kinematic bottom surface boundary condition is revelead!

$$
\begin{aligned}
&w(z) = -\mu^2 \nabla_H \cdot \int_{-h}^{z} \mathbf{u} \, dz'
- \underbrace{\left[\mu^2 \mathbf{u}(-h)(\nabla_H \cdot h) - w(-h)\right]}_{KBBC = 0}
\end{aligned}
$$

Now, we've reduced our expression for $w(z)$ to a single term. Cleaning up a bit, we will move forward with the following term for $w$:
<div class="equation-block">
  <strong>Expression for Vertical Velocity $w$</strong><br>
  $$
\begin{aligned}
&w = -\mu^2 \nabla_H \cdot \int_{-h}^{z} \mathbf{u} \, dz'
\end{aligned}
  $$
</div>

Why is this useful? Well we achieved one goal of finding some expression for $w$ that only relies on $\mathbf{u}$. And now we can substitute this into Eq. X to find the Taylor Series expansion of the horizontal velocities.

## Taylor Series Expansions of Velocities

Returning back to (), we can now write:

$$
\begin{aligned}
\frac{\partial \mathbf{u}}{\partial z} = \nabla_H w = \nabla_H\left[-\mu^2 \nabla_H \cdot \int_{-h}^{z} \mathbf{u} \, dz' \right]\\
\end{aligned}
$$

So we can substitute this into the Taylor Series expansion:

$$
\begin{aligned}
\mathbf{u} = \mathbf{u}_B + (z+h) \left(\nabla_H\left[-\mu^2 \nabla_H \cdot \int_{-h}^{z} \mathbf{u} \, dz' \right]\right) + ...\\
\end{aligned}
$$

After condensing a bit, we see that the series is already $\mathcal{O}(\mu^2)$ by the linear term.

$$
\begin{aligned}
\mathbf{u} = \mathbf{u}_B - \mu^2(z+h) \left(\nabla_H\left[ \nabla_H \cdot \int_{-h}^{z} \mathbf{u} \, dz' \right]\right) + ...\\
\end{aligned}
$$

So only including these first two terms is acceptable. Now, this term as-is is a bit useless, since $\mathbf{u}$ appears on both sides so we're defining it in terms of itself. Some final manipulations are needed to get neat functions for $\mathbf{u}$ and $w$. We begin with $w$ by substituting this into:

$$
\begin{aligned}
&w = -\mu^2 \nabla_H \cdot \int_{-h}^{z} \mathbf{u} \, dz' = -\mu^2 \nabla_H \cdot \int_{-h}^{z} \left[  \mathbf{u}_B - \mu^2(z+h) \left(\nabla_H\left[ \nabla_H \cdot \int_{-h}^{z} \mathbf{u} \, dz' \right]\right) \right] \, dz'
\end{aligned}
$$

Again this may seem circular to define $w$ in terms of $\mathbf{u}$, but if we distribute through the $\mu^2$ we see that the term containing $\mathbf{u}$ is $\mathcal{O}(\mu^4)$ and can be neglected:

$$
\begin{aligned}
&w = -\mu^2 \nabla_H \cdot \int_{-h}^{z} \mathbf{u} \, dz' = - \nabla_H \cdot \int_{-h}^{z} \left[  \mu^2\mathbf{u}_B - \underbrace{\mu^4(z+h) \left(\nabla_H\left[ \nabla_H \cdot \int_{-h}^{z} \mathbf{u} \, dz' \right]\right)}_{\mathcal{O}(\mu^4)} \right] \, dz'
\end{aligned}
$$

So we're left with the much simpler:

$$
\begin{aligned}
&w =  - \mu^2 \nabla_H \cdot \int_{-h}^{z} \mathbf{u}_B  \, dz' + \mathcal{O}(\mu^4)
\end{aligned}
$$

Finally, we're just left with the depth integral. $\mathbf{u}_B$ is *not* a function of $z$- it's just the point we're expanding about. So the integrand is constant in the dimension of integration and this is simply linear. This leads to our final expression for $w$:

$$
\begin{aligned}
&w =  - \mu^2 \nabla_H \cdot \left[\mathbf{u}_B(z+h)\right]  + \mathcal{O}(\mu^4)
\end{aligned}
$$

<div class="equation-block">
  <strong>Expression for Vertical Velocity $w$</strong><br>
$$
\begin{aligned}
&w =  - \mu^2 \nabla_H \cdot \left[\mathbf{u}_B(z+h)\right]  + \mathcal{O}(\mu^4)
\end{aligned}
$$
</div>

Finally, $\mathbf{u}$ can be found. Depth-integration is again useful. We can get this fairly directly through depth integrating $\frac{\partial\mathbf{u}}{\partial z}$ and using our expression for $w$:

$$
\begin{aligned}
& \frac{\partial\mathbf{u}}{\partial z} = \nabla_H w = \nabla_H \left[ - \mu^2 \nabla_H \cdot \left[\mathbf{u}_B(z+h)\right] \right] + \mathcal{O}(\mu^4) \\
& \int_{-h}^{z} \frac{\partial\mathbf{u}}{\partial z'} \, dz'  = \int_{-h}^{z} \nabla_H \left[ - \mu^2 \nabla_H \cdot \left[\mathbf{u}_B(z'+h)\right] \right] dz'+ \mathcal{O}(\mu^4)
\end{aligned}
$$

The LHS integral simplies nicely to an integrand of just $d\mathbf{u}$. The RHS is a bit more complicated, but still very little actually depends on $z'$. It's helpful to expand:

$$
\begin{aligned}
& \mathbf{u}(z) - \mathbf{u}(-h)  = - \mu^2 \underbrace{\int_{-h}^{z} \nabla_H \left[  \nabla_H \cdot (\mathbf{u}_Bz' ) \right]\, dz'}_{(1)} - \mu^2 \underbrace{\int_{-h}^{z} \nabla_H \left[\nabla_H \cdot (\mathbf{u}_Bh) \right]dz'}_{(2)}+ \mathcal{O}(\mu^4)
\end{aligned}
$$

Integral (2) has no $z$-dependency in the integrand at all, so the result is simply linear. Here we'll also swap out notation on the LHS

$$
\begin{aligned}
& \mathbf{u} - \mathbf{u}_B  = - \mu^2 \underbrace{\int_{-h}^{z} \nabla_H \left[  \nabla_H \cdot (\mathbf{u}_Bz' ) \right]\, dz'}_{(1)} - \mu^2(z+h) \nabla_H \left[\nabla_H \cdot (\mathbf{u}_Bh) \right]dz'+ \mathcal{O}(\mu^4)
\end{aligned}
$$

Integral (1) is a bit trickier, but note first how $z'$ can be brought outside the $\nabla_H$ operator (as either a divergence or gradient) since $\nabla_H$ only acts in the horizontal:

$$
\begin{aligned}
& \int_{-h}^{z} \nabla_H \left[  z' \nabla_H \cdot \mathbf{u}_B \right]\, dz' = \int_{-h}^{z} z' \nabla_H \left[   \nabla_H \cdot \mathbf{u}_B \right]\, dz'
\end{aligned}
$$

And now everything to the right of the first gradient doesn't depend on $z'$, so it can be taken out and the integral is simply evaluated:

$$
\begin{aligned}
&  \int_{-h}^{z} z' \nabla_H \left[   \nabla_H \cdot \mathbf{u}_B \right]\, dz' = \nabla_H \left[   \nabla_H \cdot \mathbf{u}_B \right]\int_{-h}^{z} z' \, dz' = \nabla_H \left[   \nabla_H \cdot \mathbf{u}_B \right] \left(\frac{z^2}{2}-\frac{h^2}{2}\right)
\end{aligned}
$$

Substituting this into our result from earlier we have:

$$
\begin{aligned}
& \mathbf{u} - \mathbf{u}_B  = - \mu^2 \left[\nabla_H \left[   \nabla_H \cdot \mathbf{u}_B \right] \left(\frac{z^2}{2}-\frac{h^2}{2}\right)\right]- \mu^2(z+h) \nabla_H \left[\nabla_H \cdot (\mathbf{u}_Bh) \right]+ \mathcal{O}(\mu^4)
\end{aligned}
$$

And rearranging a bit gives us our expression for $\mathbf{u}$:


<div class="equation-block">
  <strong>Expression for Horizontal Velocity $\mathbf{u}$</strong><br>
  $$
\begin{aligned}
&\mathbf{u} = \mathbf{u}_B  + \mu^2 \left[\left(\frac{h^2}{2}-\frac{z^2}{2}\right)\nabla_H \left[   \nabla_H \cdot \mathbf{u}_B \right] \right]  - \mu^2(z+h) \nabla_H \left[\nabla_H \cdot (h\mathbf{u}_B) \right]+ \mathcal{O}(\mu^4)
\end{aligned}
  $$
</div>


## Extension to an Arbitrary Depth
Starting with Nwogu (1993), it is noted that the reference elevation need not be $z_\alpha=-h$ and the extension to these arbitrary elevations is straightforward. In fact, most
FUNWAVE-like models choose a $z_\alpha$ to optimize dispersion behavior or some other desired model
quality.

To see this, let's start by taking a Taylor Series expansion about the bed up to $z_\alpha$. (To keep things concise, we'll just use the standard Taylor Series form, assuming we don't yet know the derivative)

$$
\begin{aligned}
\mathbf{u}_\alpha = \mathbf{u}_B + \underbrace{(z_\alpha + h)\left(\frac{\partial \mathbf{u}}{\partial z}\right)_B}_{\mathcal{O}(\mu^2)} + \mathcal{O}(\mu^4)
\end{aligned}
$$

The second term on the RHS is $\mathcal{O}(\mu^2)$, and we ultimately neglect the $\mathcal{O}(\mu^4)$ terms. So $\mathbf{u}_\alpha$ and $\mathbf{u}_B$ are equivalent up to second order:

$$
\begin{aligned}
\mathbf{u}_B = \mathbf{u}_\alpha - \mathcal{O}(\mu^2)
\end{aligned}
$$

If we substitute this expression for $\mathbf{u}_B$ into our velocity expansions about the bed, we see that we quickly accumulate  $\mathcal{O}(\mu^4)$ terms as we bring the $\mu^2$ in. We begin with $w$

$$
\begin{aligned}
&w =  - \mu^2 \nabla_H \cdot \left[\mathbf{u}_B(z+h)\right]  + \mathcal{O}(\mu^4) =  - \mu^2 \nabla_H \cdot \left[\left( \mathbf{u}_\alpha - \mathcal{O}(\mu^2)\right)(z+h)\right]  + \mathcal{O}(\mu^4) \\
\end{aligned}
$$

If we distribute in we see that we get a $\mathcal{O}(\mu^4)$ term

$$
\begin{aligned}
&w =   -  \nabla_H \cdot \left[\left( \mu^2\mathbf{u}_\alpha - \mathcal{O}(\mu^4)\right)(z+h)  \right]  + \mathcal{O}(\mu^4) \\
& w =   -\mu^2 \nabla_H \cdot \left[\mathbf{u}_\alpha(z+h)\right]  + \mathcal{O}(\mu^4) \\
\end{aligned}
$$

So for $w$, the expression is the same whether we use  or $\mathbf{u}_\alpha$ or $\mathbf{u}_B$. Since we derived $\mathbf{u}$ from $\frac{\partial u}{\partial z}=\nabla_H w$, we can substitute velocities there as well. Thus, in summary, we have:



<div class="equation-block">
  <strong>Taylor Series Velocity Expansions about $z_\alpha$</strong><br>
  $$
\begin{aligned}
&\mathbf{u} = \mathbf{u}_\alpha  + \mu^2 \left[\left(\frac{h^2}{2}-\frac{z^2}{2}\right)\nabla_H \left[   \nabla_H \cdot \mathbf{u}_\alpha \right] \right]  - \mu^2(z+h) \nabla_H \left[\nabla_H \cdot (h\mathbf{u}_\alpha) \right]+ \mathcal{O}(\mu^4) \\
& w =   -\mu^2 \nabla_H \cdot \left[\mathbf{u}_\alpha(z+h)\right]  + \mathcal{O}(\mu^4)
\end{aligned}
  $$
</div>


## Implications of the Taylor Series Expansion
With these two expansions completed, it is worth considering what they imply about the actual solution. As a Boussinesq model, FUNWAVE-TVD truncates the vertical structure via a Taylor Series to a given order, in this case, $\mu^2$. By looking at 

## Conclusion
This concludes the Taylor Series expansions of the velocity terms, leaving only pressure remaining. Pressure ends up being substantially more complicated. This is perhaps unsurprising, since non-hydrostatic models often spend much of their time calculating pressures. 

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
    <div class="forward-link"><a href="/notes/FWguide/pressure/">Taylor Series Expansion of Pressure</a></div>
  </div>
</div>
