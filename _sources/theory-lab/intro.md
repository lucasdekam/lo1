# Introduction

In this practical we will look at the oxygen reduction reaction. This half-reaction occurs at the cathode of a fuel cell. Fuel cells are used to convert hydrogen gas and oxygen gas into water and energy in hydrogen-powered cars like the Toyota Mirai, for example. However, the oxygen reduction reaction is also interesting as a simple reaction that teaches us a lot about reactions at electrodes in general.

In the oxygen reduction reaction, an oxygen molecule gains four electrons and four protons to form water:

$$\mathrm{O_2 + 4 H^+ + 4e^- \to 2 H_2O}.$$

The electrons come from the electrode. The protons come from the electrolyte (we consider an acidic electrolyte here). 

The transfer of protons and electrons does not happen all at once; the reaction consists of multiple steps. Although multiple reaction mechanisms are possible, we will look at the simple mechanism from the work of {cite:t}`norskov2004origin`. In the first step, an oxygen molecule adsorbs on the surface and dissociates:

$$\mathrm{O_2 + 2 * \to 2 O^*}.$$

Here $*$ is an empty adsorption site on the surface and $\mathrm{O^*}$ is an adsorbed oxygen atom. Next, a proton and an electron is transferred to each $\mathrm{O^*}$:

$$\mathrm{O^* + H^+ + e^- \to HO^*}$$

and this can happen one more time to form water:

$$\mathrm{HO^* + H^+ + e^- \to H_2O + *}.$$

The electrode surface acts as a catalyst: it binds the *intermediates* $\mathrm{O^*}$ and $\mathrm{HO^*}$ during the reaction. The reduction of oxygen to water causes a flow of electrons: an electric current. The field of science that studies catalytic reactions that drive or require an electric current is called *electrocatalysis*. 

## The binding energy in electrocatalysis

In the lab part, you will find or have found that different catalysts (different electrodes) have a different activity for the oxygen reduction reaction. 


Gibbs free energy of adsorption and computational hydrogen electrode.



## Quantum chemistry

The following text is inspired by the lecture notes of Christopher Stein (TU Munich) and Katharina Doblhoff-Dier (Leiden Univ.), and the Wikipedia page on [density functional theory](https://en.wikipedia.org/wiki/Density_functional_theory). 

In this practical course we consider electrocatalysis at the scale of atoms. In particular, we are interested in the chemical bond between a reaction intermediate and a catalyst. A chemical bond forms due to the interaction between the electrons of the reaction intermediate and the catalyst. So: to understand chemical bonds, we need to understand electrons.

Electrons are very small and very light, much lighter even than atomic nuclei. Electrons do not obey Newton's laws of classical mechanics. Instead, we need to consider quantum mechanics. The central quantity in quantum mechanics is the wave function. If we know the wave function of something, we can use it to calculate any observable quantity.

The wave function of all electrons in a certain system (for example, a molecule) can be found by solving the Schrödinger equation,

$$ \hat{H}(\mathbf{R}) \Psi(\mathbf{R}) = E \Psi(\mathbf{R}) $$

where $\mathbf{R}$ denotes all coordinates of all electrons in the system, $\hat{H}$ is the Hamiltonian of the system, $E$ is the energy and $\Psi$ is the wavefunction. Mathematically, $\Psi(\mathbf{R})$ is an *eigenfunction* of $\hat{H}$, and $E$ the corresponding *eigenvalue*. 

```{note} 
:class: dropdown
This is the time-independent Schrödinger equation. We don't need to consider dynamics of electrons. We also assume that all nuclei are fixed in space. We can make this assumption because the nuclei are much heavier than the electrons. This approximation is also called the Born-Oppenheimer approximation.
```

The Hamiltonian is an operator. You can think of it as a 'recipe' that tells you how the wave function is related to the energy. For the electrons in a molecule, the Hamiltonian is the sum of several parts:

$$ \hat{H} = \hat{T} + \hat{V}_{ne} + \hat{V}_{ee}.$$

The first term, $\hat{T}$, is the kinetic energy of the electrons; $\hat{V}_{ne}$ the interaction between the electrons and the nuclei; and $\hat{V}_{ee}$ is the interaction between the different electrons. The kinetic energy is no problem; there is a closed-form formula ($\hat{T} = -\hbar^2/2m \times  \nabla^2 \Psi$ for each electron). The interaction between electrons and nuclei is also no problem, because we assumed that the nuclei are fixed in space. However, the electron-electron interaction is very complicated (read about the 'many-body problem' on [Wikipedia](https://en.wikipedia.org/wiki/Many-body_problem)). 

To solve the Schrödinger equation for many-electron systems (where many is more than three) on a computer within a reasonable calculation time, people have come up with many different approximations, from [Hartree-Fock theory](https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method) to [coupled-cluster theory](https://en.wikipedia.org/wiki/Coupled_cluster). All these methods revolve around calculating the wave function. Some of these methods can be extremely accurate, but it takes a lot of time to do such calculations on a computer. It turns out that there is another approach that achieves decent accuracy in much less computation time: the *density functional theory* (DFT). 

In density functional theory, the central quantity is the electron density, $\rho(\mathbf{r})$, where $\mathbf{r}$ is the vector of spatial coordinates. Walter Kohn and Pierre Hohenberg showed that $\rho(\mathbf{r})$ also contains all relevant information about the system. Later, Walter Kohn and Lu Jeu Sham showed that the problem of calculating the electronic density can be reduced to solving the Schrödinger equation of noninteracting particles (so no $\hat{V}_{ee}$), which made the method very easy to use, while still being accurate. Walter Kohn was awarded the [Nobel Prize in 1998](https://www.nobelprize.org/prizes/chemistry/1998/summary/). 

The method of Kohn and Sham comes down to solving $N$ Schrödinger equations, one for each electron. The electrons move in an effective potential $V_\mathrm{eff}(\mathbf{r})$, which depends on the electron density and the positions of the nuclei. The system of these particular $N$ Schrödinger equations is known as the Kohn-Sham equations; they read

$$ \left( -\frac{\hbar^2}{2m} \nabla^2 + V_\mathrm{eff}(\mathbf{r}) \right) \psi_i(\mathbf{r}) = \varepsilon_i \psi_i(\mathbf{r}). $$

Each electron has its own one-electron wavefunction $\psi_i$, known as an *orbital*. These orbitals have an associated energy $\varepsilon_i$. The density is calculated from the orbitals as

$$ \rho(\mathbf{r}) = \sum_{i=1}^{N} |\psi_i(\mathbf{r})|^2. $$

The effective potential depends on the density as follows:

$$ V_\mathrm{eff}(\mathbf{r}) = V_\mathrm{ext}(\mathbf{r}) + \int \frac{\rho(\mathbf{r}')}{4\pi\varepsilon_0 |\mathbf{r} - \mathbf{r}'|} \mathrm d \mathbf{r}' + \frac{\delta E_\mathrm{xc}[\rho(\mathbf{r})]}{\delta \rho(\mathbf{r})}.$$

For this course, it is not important to understand what all these terms mean in detail. In simple words:

* $V_\mathrm{ext}$ (ext for external) describes the electrostatic interaction of the electron with the various nuclei in the system. 

* The integral of $\rho$ describes the electrostatic interaction of the electron with the negatively charged electron cloud, consisting of all electrons in the system. 

* The *exchange-correlation functional* $E_\mathrm{xc}$ contains corrections to the approximations made by Kohn and Sham. Finding better exchange-correlation functionals is an ongoing topic in research. Different systems might be described better by different exchange-correlation functionals.

```{note} 
:class: dropdown
In case you are wondering what the symbols in the last term mean, here's an explanation. A *function* $f(x)$ takes a number as input, and outputs a number. For example: $f(x)=x^2$ turns $x=2$ into $f(2)=4$. A *functional* $F[f(x)]$ takes an entire function as input, and outputs a number. Functionals also have their own kind of derivative, denoted with a $\delta$. 
```

It follows that the Kohn-Sham equations depend on $\rho$ through $V_\mathrm{eff}$. To write down and solve the Kohn-Sham equations, we therefore need to know $\rho$, but to know $\rho$ we need to know the $\psi_i$, which we get by solving the Kohn-Sham equations. How do we solve this circular problem? Luckily, it turns out that when we solve the Kohn-Sham equations for a certain $\rho$, we get out $\psi_i$ that give a 'better' $\rho$. 'Better' means: lower energy, because we are looking for the lowest energy state (ground state). This state gives us the most optimal configuration. By starting out with a certain guess for $\rho$, we can keep solving the Kohn-Sham equations with an improved $\rho$ every time. When our energy does not change anymore, we know that we are at a minimum: we are *converged*. Such a 'circular' procedure is called a *self-consistent procedure*. 

Now the question remains: how do we solve the Kohn-Sham equations numerically, and what are the $\psi_i$? Depending on whether you have to simulate metal electrodes or porphyrin molecules, you will use different representations of $\psi_i$. These representations will be discussed in the next pages, and you will use software to use the Kohn-Sham equations self-consistently yourself.