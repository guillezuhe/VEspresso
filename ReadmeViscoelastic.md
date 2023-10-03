# Viscoelastic Espresso

This version of espresso introduces a non-newtonian carrier medium using the non-linear Generalized Langevin Equation.

## Generalized Langevin Equation

Following *Andrew D. Baczewski1 and Stephen D. Bond Numerical Integration of the Extended Variable Generalized Langevin Equation with a
Positive Prony Representable Memory Kernel*.

The **General Langevin Equation (GLE)** explains the movement of a particle for more general viscoelastic fluids and solids. We thus assume no specific form of the resistance or memory function $\Gamma (t-t')$, other than causality.

$$
\begin{cases}
m \frac{d\bold{v}}{dt} = \bold{F^c}(\bold{x}(t)) + \bold{F^r}(t) - \int_{0}^{t}{\Gamma (t-t') \bold{v}(t')dt'} \\
\bold{v(t)} = \frac{d\bold{x}(t)}{dt}
\end{cases}
$$

with initial conditions $\bold{x}(0) = \bold{x_0}$ and $\bold{v}(0) = \bold{v_0}$. $F^c$ is a conservative force and $F^r$ is a random force, that must satisfy the Fluctuation-Dissipation Theorem (FDT):

$$
\langle F_i^r(t') F_j^r(t) \rangle = k_B T \; \Gamma (t-t') \; \delta_{ij}
$$

Subscripts $i$ and $j$ denote the different components of the vector.

In order to solve this general equation with correlated noise, we assume that the memory kernel function can be represented as a Prony series of decaying exponentials with $N_k$ terms.

$$
\Gamma (t) = \sum_{k=1}^{N_k}{\frac{c_k}{\tau_k} \exp[-t / \tau_k]} \; \; \; \;\text{ for } t \ge 0
$$

where $c_k$ and $\tau_k$ are fitting parameters. We will discuss in following sections how to obtain these fitting parameters. We can now define an extended variable $Z_{i,k} (t)$:

$$
Z_{i,k}(t) = - \int_{0}^{t}{\frac{c_k}{\tau_k} \exp[- (t - t') / \tau_k] v_i (t') dt'} 
$$

and making use of $Z_{i,k} (t)$ we can rewrite the GLE for each spatial component:

$$
\begin{cases}
m dv_i = F_i^c (\bold{x}(t))dt + \sum_{k=1}^{N_k}{Z_{i,k} (t) dt} + F_i^r (t) dt \\
dx_i (t) = v_i (t) dt
\end{cases}
$$

Now, instead of writing $Z_{i,k} (t)$ in integral form, we can write it in differential form to obtain a Stochastical Differential Equation (SDE).

$$
dZ_{i,k} (t) = - \frac{1}{\tau_k} Z_{i,k} (t) dt - \frac{c_k}{\tau_k} v_i (t) dt
$$

This way, we can couple this equation with the previous set and solve them as a system, obtaining $x_i$, $v_i$ and $z_{i,k}$ without the need of computing the integral.

Let us now face the contribution of the random force $F^r$ to make it obey the FDT. To this end, we consider the following SDE

$$
dF_{i,k} (t) = - \frac{1}{\tau_k} F_{i,k} (t) dt + \frac{1}{\tau_k} \sqrt{2 k_B T c_k} dW_{i,k} (t)
$$

If $W_{i,k}$ is a standard Wiener process, this SDE corresponds to an Ornstein-Uhlenbeck (OU) process. Using the OU properties, one can see that:

$$
\langle F_{i,k} (t') F_{i,k} (t) \rangle = k_B T \frac{c_k}{\tau_k} \exp[-(t-t') / \tau_k]
$$

and we can straightforwardly write:

$$
F_i^r (t) = \sum_{k=1}^{N_k}{F_{i,k} (t)}
$$

Combining both results, we can define the final extended variable $S_{i,k} (t) = Z_{i,k} (t) + F_{i,k} (t)$ and achieve a final expression for our set of equations:

$$
\begin{cases}
m dv_i = F_i^c (\bold{x}(t))dt + \sum_{k=1}^{N_k}{S_{i,k} (t) dt} \\
dx_i (t) = v_i (t) dt \\
dS_{i,k} (t) = - \frac{1}{\tau_k} S_{i,k} (t) dt - \frac{c_k}{\tau_k} v_i (t) dt + \frac{1}{\tau_k} \sqrt{2 k_B T c_k} dW_{i,k} (t)
\end{cases}
$$

**Note**: It can be shown, that in the limit of small $\tau_k$, the GLE converges to the Newtonian Langevin equation (see original paper).