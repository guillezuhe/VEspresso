# Viscoelastic Espresso

This version of espresso introduces a non-newtonian carrier medium using the non-linear Generalized Langevin Equation.

## Generalized Langevin Equation

Following *Andrew D. Baczewski1 and Stephen D. Bond Numerical Integration of the Extended Variable Generalized Langevin Equation with a
Positive Prony Representable Memory Kernel*.

The **General Langevin Equation (GLE)** explains the movement of a particle for more general viscoelastic fluids and solids. We thus assume an arbitrary form of the resistance which will be given by a typical Newtonian component plus a viscoelastic contribution determined by its memory function $\Gamma (t-t')$.

$$
\begin{cases}
m \frac{d\bold{v}}{dt} = \bold{F^c}(\bold{x}(t)) -\zeta \bold{v}(t) + \bold{F^r_N} (t) - \int_{0}^{t}{\Gamma_m (t-t') \bold{v}(t')dt'} + \bold{F^r}(t) \\
\bold{v(t)} = \frac{d\bold{x}(t)}{dt}
\end{cases}
$$

with initial conditions $\bold{x}(0) = \bold{x_0}$ and $\bold{v}(0) = \bold{v_0}$. $F^c$ is a conservative force and. $F^r_N$ and $F^r$ are Newtonian and viscoelastic random forces respectively, that must satisfy the Fluctuation-Dissipation Theorem (FDT):

$$
\langle F_{N, i}^r(t') F_{N, j}^r(t) \rangle = 2 \zeta k_B T \; \delta(t-t') \; \delta_{ij}
$$

$$
\langle F_i^r(t') F_j^r(t) \rangle = k_B T \; \Gamma_m (t-t') \; \delta_{ij}
$$

Subscripts $i$ and $j$ denote the different components of the vector.

In order to solve this general equation with correlated noise, we assume that the memory kernel function can be represented as a Prony series of decaying exponentials with $N_k$ terms.

$$
\Gamma_m (t) = \sum_{m=1}^{N_m}{\Kappa_m \exp[-t / \tau_m]} \; \; \; \;\text{ for } t \ge 0
$$

where $c_k$ and $\tau_k$ are fitting parameters. We will discuss in following sections how to obtain these fitting parameters. We can now define an extended variable $Z_{i,k} (t)$:

$$
Z_{i,m}(t) = - \int_{0}^{t}{\exp[- (t - t') / \tau_m] v_i (t') dt'} 
$$

and making use of $Z_{i,k} (t)$ we can rewrite the GLE for each spatial component:

$$
\begin{cases}
m dv_i = F_i^c (\bold{x}(t))dt -\zeta v_i (t) dt + F^r_{N,i} (t) dt + \sum_{m=1}^{N_m}{\Kappa_m Z_{i,m} (t) dt} + F_i^r (t) dt \\
dx_i (t) = v_i (t) dt
\end{cases}
$$

Now, instead of writing $Z_{i,k} (t)$ in integral form, we can write it in differential form to obtain a Stochastical Differential Equation (SDE).

$$
dZ_{i,m} (t) = - \frac{1}{\tau_m} Z_{i,m} (t) dt - v_i (t) dt
$$

This way, we can couple this equation with the previous set and solve them as a system, obtaining $x_i$, $v_i$ and $z_{i,k}$ without the need of computing the integral.

Let us now face the contribution of the random force $F^r$ to make it obey the FDT. To this end, we consider the following SDE

$$
df_{i,m} (t) = - \frac{1}{\tau_m} f_{i,k} (t) dt + \sqrt{2 k_B T / \zeta_m} dW_{i,m} (t)
$$

If $W_{i,m}$ is a standard Wiener process, this SDE corresponds to an Ornstein-Uhlenbeck (OU) process. Using the OU properties, one can see that:

$$
\langle f_{i,m} (t') f_{i,m} (t) \rangle = \frac{k_B T}{\Kappa_m} \exp[-(t-t') / \tau_m]
$$

and we can straightforwardly write:

$$
F_i^r (t) = \sum_{k=m}^{N_m}{\Kappa_m f_{i,m} (t)}
$$

Combining both results, we can define the final extended variable $U_{i,m} (t) = Z_{i,m} (t) + f_{i,m} (t)$ and achieve a final expression for our set of equations:

$$
\begin{cases}
m dv_i = F_i^c (\bold{x}(t))dt -\zeta v_i (t) dt + F^r_{N,i} (t) dt + \sum_{m=1}^{N_m}{\Kappa_m U_{i,m} (t) dt} \\
dx_i (t) = v_i (t) dt \\
dU_{i,m} (t) = - \frac{1}{\tau_m} U_{i,m} (t) dt - v_i (t) dt + \sqrt{2 k_B T / \zeta_m} dW_{i,k} (t)
\end{cases}
$$

**Note**: It can be shown, that in the limit of small $\tau_k$, the GLE converges to a typical Langevin equation (see original paper).

### Physical interpretation

The simplest way of modelling the viscoelastic behavior of non-Newtonian fluids is their treatment as combinations of two basic elements: 

* Dampers, representing viscous friction, which exert a force $f_{\zeta} = - \zeta \dot{x}$.
  
* Springs, standing for elastic restoring forces, with a general form $f_{\Kappa} = - \Kappa x$


As seen in the previous derivation, our memory function must be decomposed in a series of decaying exponentials. A single exponential kernel function, without the extra Newtonian contribution, is one of the simplest viscoelastic models, a Maxwell fluid, which corresponds to a viscous damper and a elastic element connected sequentially.

$$
\Gamma (t) = \frac{\zeta_m}{\tau_m} \exp{(-t / \tau_m)} 
$$

where $\zeta_m$ is the Maxwellian viscous friction coefficient, and $\tau_m$ is the relaxation time, the ratio of viscous and elastic coefficients $\tau_m = \zeta_m / \Kappa_m$.

The addition of the Newtonian element leads to the so called Jeffreys model: a Maxwell chain connected in parallel with a purely viscous element.

$$
\Gamma (t) = 2 \, \zeta \, \delta (t) + \frac{\zeta_m}{\tau_m} \exp{(-t / \tau_m)} 
$$

This way, expressing our memory function as sum of exponentials takes us to a General Jeffreys model, that is, a Newtonian viscous element to an arbitrary number $N_m$ of Maxwell elements.

$$
\Gamma(t) = 2 \, \zeta \, \delta (t) + \sum_{m=1}^{N_m}{\frac{\zeta_m}{\tau_m} \exp{(-t / \tau_m)}}
$$

Therefore, the Generalized Langevin Equation approach is exact when the memory kernel is a sum of exponential, but provides a very efficient and accurate method to treat arbitrary memory functions as long as they can be fitted to an exponential sum. This decomposition is specially useful, for example, in power-law scaling memory kernels.

$$
\Gamma (t) = K t^{-n}
$$

This power-law rheology is usually observed in many complex materials such as biofluids, foods, cross-linked polymers, microgels, and hydrogels.

Another typical memory kernel is the Rouse model, another bead-spring model representative of diluted flexible polymers in Newtonian solvents. This model is just a sum of exponentials over the relaxation of the chain's normal modes, without any Newtonian contribution:

$$
\Gamma = n k_B T \zeta \sum_{p=1}^{N} \exp{(-t / \tau_p)}
$$

### Equation integration

We assume a constant timestep $\Delta t$ and we will use the notation $x_i (n\Delta t) = x_i^n$. Therefore, having the values $x_i^n$, $v_i^n$ and $s_{i,k}^n$, we will update to the $(n+1)$th time step following this method:

1. Update $v_i$ by a half step:
   $$
   v_i^{n+1/2} = v_i^n + \frac{\Delta t}{2 m} F_i^c(x^n) + \frac{\Delta t}{2 m} \left( -\zeta v_i^n + \sqrt{\frac{2 k_B T \zeta}{\Delta t}} B_{i,N}^n + \sum_{m=1}^{N_m}{\Kappa_m U_{i,m}^n} \right)
   $$
2. Update $x_i$ by a full step:
   $$
   x_i^{n+1} = x_i^n + \Delta t \; v_i^{n+1/2}
   $$
3. Update the forces and the viscoelastic contribution $U_{i,m}$ a by a full step with the new positions and velocities:
   $$
   U_{i,m}^{n+1} = U_{i,m}^n - \left(\frac{1}{\tau_m} U_{i,m}^n + v_i^{n+1/2} - \sqrt{\frac{2 k_B T}{\zeta_m \, \Delta t}} B_{i,m}^n \right) \Delta t
   $$
4. Update $v_i$ by another half step:
   $$
   v_i^{n+1} = v_i^{n+1/2} + \frac{\Delta t}{2 m} F_i^c(x^{n+1}) + \frac{\Delta t}{2 m} \left( -\zeta v_i^{n+1} + \sqrt{\frac{2 k_B T \zeta}{\Delta t}} B_{i,N}^{n+1} + \sum_{m=1}^{N_m}{\Kappa_m U_{i,m}^{n+1}} \right)
   $$

$B_{i,N}^n$ and $B_{i,m}^n$ are uncorrelated random numbers chosen from a Gaussian distribution of zero mean and variance unity $N(0,1)$.


### Non-linearity

All the previously discussed viscoelastic models, are restricted within the linear response limit, but many interesting properties emerge in when considering the non-linear rheology of materials. Among these properties, we can point shear rate-dependent viscosity (shear-thinning or thickening), yield stresses or even normal stress differences.

To implement these non-linearities in the Langevin Equation, a few models have been developed. Some of them, comprise the introduction of negative memory functions [[Ref]](https://doi.org/10.1038/s41467-018-03345-2), but it is a purely phenomenologycal fitting approach. Here we will follow the approach introduced by [I. Goychunk](https://doi.org/10.1073/pnas.2205637119), which is a non-linear Generalized Jeffreys-Langeving model, where the spring and damper constants are velocity dependent: $\Kappa_m (v)$ and $\zeta_m (t)$.

$$
\begin{cases}
\bold{v(t)} = \frac{d\bold{x}(t)}{dt} \\
m \frac{d\bold{v} (t)}{dt} = \bold{F^c}(\bold{x}(t)) -\zeta \bold{v}(t) + \bold{F^r_N} (t) + \sum_{m=1}^{N_m}{\Kappa_m (v) \, \bold{U_{m}} (t)} \\
\frac{d\bold{U_{m}} (t)}{dt} = - \frac{1}{\tau_m (v)} \bold{U_{m}} (t) - \bold{v} (t) + \sqrt{2 k_B T / \zeta_m (v)} \bold{\xi_m} (t)
\end{cases}
$$

where it is shown that the relaxation time can also be velocity dependent. A basic case of this model, implemented in this version of espresso, consists on a constant relaxation time $\tau_m (v) = \tau_m$, and a friction coefficient following a particular case of the Carreau-Yasuda model.

$$
\zeta_m (v) = \frac{\zeta_{m0}}{(1 + |v/v_c|^b)^a}
$$

where the model parameters are: $\zeta_{m0}$, the zero-shear viscous friction coefficient, $v_c$, the critical velocity that marks the transition between constant and power-law viscosity, $a$ and $b$ the scaling exponents of the model.

Lastly, it can be shown that for vanishingly small relaxation times $\tau_m \rightarrow 0$, this model reduces to a classical Langevin equation with a viscosity that is non-linear in velocity:

$$
m \frac{d\bold{v} (t)}{dt} = \bold{F^c}(\bold{x}(t)) -\zeta_\text{{eff}} (v) \bold{v}(t) + \bold{F^r_N} (v,t)
$$

with a friction coefficient 

$$
\zeta_\text{{eff}} (v) = \zeta + \sum_{m=1}^{N_m} {\zeta_m (v)}
$$



## Espresso implementation

This espresso version works in the same fashion of the main version. To include carrier medium viscoelasticity, we just need to define the involved parameters and pass them to espresso through the python interface.

The viscoelastic parameters have been introduced in espresso as particle properties so that they can be individually defined. However, for monodisperse particles in isotropic mediums, parameters are expected to take the same value for all the particle set.