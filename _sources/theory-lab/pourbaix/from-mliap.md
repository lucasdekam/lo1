# Pourbaix diagram from calculations

The final task is to use the energies we calculated from MACE to make a Pourbaix diagram. We need to connect energies from an atomistic calculation to the real world. This connection is provided by thermodynamics. The method here is mostly based on {cite:t}`persson2012prediction`.

Check out the Materials Project Pourbaix diagram generator at https://next-gen.materialsproject.org/pourbaix. It generates Pourbaix diagrams based on a database of PBE-DFT calculations. The challenge is: can we do better than these DFT calculations (which required a supercomputer) with just our laptop and the MACE foundational model?

## From electronic energy to chemical potential

From DFT, or MACE, we get an electronic energy $E^\mathrm{DFT}$. We can also compute the energy due to vibration of the nuclei using the forces on the nuclei we get from DFT calculations, but this contribution is usually small. Hence, we approximate the internal energy as

$$U \approx E^\mathrm{DFT}.$$

We saw before that we can then calculate the Gibbs energy $G=U-TS+pV$. The $pV$ contributions are usually very small, so we approximate these to be zero. The chemical potential $\mu$ is the Gibbs energy per particle (for molecules) or per structural unit (in a periodic lattice). We can calculate it as

$$\mu \approx E^\mathrm{DFT} - Ts$$

where $s$ is the entropy per particle and $T$ the temperature in Kelvin.

But we have a problem: the chemical potential of an element in its standard state is zero. But $E^\mathrm{DFT}$ is not zero! To solve this, we introduce **reference chemical potentials**. The reference state of an element is its stable state at standard conditions. For example, a bcc structure for an iron crystal, and gas for an oxygen molecule.

When calculating the chemical potential of a material/molecule, we always **subtract the reference chemical potentials** for all the elements in the material/molecule. Sounds confusing, so let's see some examples.

## Solid elements

For all solid metal elements (pure iron, pure cobalt, ...) that we encounter in this lab practical, we are interested in their stable state at standard conditions. This state is the same as the reference state. Also, for all solids, we approximate the entropy $s \approx 0$.

The reference chemical potential is

$$\mu_\mathrm{M}^\mathrm{ref}= E^\mathrm{DFT}_\mathrm{M}-Ts_i \approx E^\mathrm{DFT}_\mathrm{M} $$

The chemical potential at standard conditions is then

$$\mu_\mathrm{M}^0 = E^\mathrm{DFT}_\mathrm{M}-Ts_\mathrm{M}^0 - \mu_\mathrm{M}^\text{ref} \approx E^\mathrm{DFT}_\mathrm{M} - E^\mathrm{DFT}_\mathrm{M} = 0$$

By our definition, the standard chemical potential of a solid element at standard conditions is indeed zero.

## Oxygen gas

At standard conditions, oxygen is a gas, so its reference chemical potential is expressed using the gaseous state. Since we calculate the energy of $\mathrm{O_2}$, we need to divide $E^\mathrm{DFT}_\mathrm{O_2}$ by $2$ to get the energy of a single oxygen atom at standard state.

It's possible to calculate $s_\mathrm{O}$ with DFT or MACE, but we just use the experimental value from {cite:t}`persson2012prediction`: $s_\mathrm{O}^0=+10.6 \times 10^{-4} \mathrm{eV}$ per oxygen atom.

The reference chemical potential is

$$\mu_\text{O}^\text{ref}\approx\frac12 E_\mathrm{O_2}^\text{DFT} - Ts_\text{O}^\text{0,exp}$$

As oxygen is an element, $\mu_\mathrm{O}^0=\frac12 E_\mathrm{O_2}^\text{DFT} - Ts_\text{O}^\text{0,exp}-\mu_\mathrm{O}^\mathrm{ref}=0.$

<!-- Common DFT methods are pretty good at calculating differences between solids, but not so good at calculating energies of molecules. So to get a more accurate Pourbaix diagram, we will make some corrections. We define the reference state of oxygen as

$$\mu_\text{O}^\text{ref}=E_\text{O}^\text{DFT}+ \Delta E^\text{correction}_\text{O} - Ts_\text{O}^\text{0,exp} + pV$$

Persson et al. use $\Delta E_\text{O}^\text{correction}=1.36$ eV/O, and $Ts_\text{O}^\text{0,exp}=0.317$ eV/O. However, we can also use $\Delta E_\text{O}^\text{correction}$ as a fitting parameter for our results to experimental formation energies. -->


## Solid oxide compounds

Up until now, all standard chemical potentials were just zero, but things get more interesting when we mix elements. For a metal oxide,

$$\mu_\mathrm{M_{x} O_{y}}^0 = E_\mathrm{M_{x} O_{y}}^\mathrm{DFT}- x\mu_\mathrm{M}^\mathrm{ref}-y\mu^\mathrm{ref}_\mathrm{O}$$

Note that we subtract the reference chemical potentials $\mu_i^\mathrm{ref}$, not standard chemical potentials $\mu_i^0$. Otherwise nothing would happen, because $\mu_i^0=0$ for elements.

We just calculated the formation energy of a metal oxide from quantum mechanics!

`````{admonition} Task
:class: tip

Calculate $\mu^0$ for your metal oxides from the MACE energies. Make a plot with $\mu^0$ from MACE on the $x$-axis, and $\mu^0$ from literature (which you gathered in the first part) on the $y$-axis. Do the data points lie on the line $x=y$? What does this mean?
`````

## Water

We use the experimental standard formation energy of water,

$$\mu^0_\mathrm{H_2O}=-2.46 \;\mathrm{eV/H_2O}$$

This can be calculated from [water data](https://en.wikipedia.org/wiki/Water_(data_page)), [hydrogen data](https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=1) and [oxygen data](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Units=SI&Mask=1#Thermo-Gas).

## Hydrogen gas

We could also have calculated the chemical potential of water with DFT and reference it to the chemical potentials of the elements:

$$\mu_\mathrm{H_2O}^0=E_\mathrm{H_2O}^\mathrm{DFT}-Ts^0_\mathrm{H_2O}+pV-2\mu_\mathrm{H}^\mathrm{ref}-\mu_\mathrm{O}^\mathrm{ref}.$$

In this equation, $\mu_\mathrm{H}^\mathrm{ref}$ is the only unknown, because we defined the chemical potential of water as the experimental value. We can solve for it as

$$
\mu_\mathrm{H}^\mathrm{ref} = \frac12 \left[ E_\mathrm{H_2O}^\mathrm{DFT}-Ts^0_\mathrm{H_2O}-\mu_\mathrm{O}^\mathrm{ref}
-\mu^0_\mathrm{H_2O} \right].
$$

You'll need the entropy of water: $s_\mathrm{H_2O}^0=+7.25 \times 10^{-4} \mathrm{eV}$ per water molecule (see [Water data](https://en.wikipedia.org/wiki/Water_(data_page))).

Of course, as it is an element, $\mu_\mathrm{H}^0=0$.

`````{admonition} Task
:class: tip

Use the reference potential for hydrogen to calculate $\mu^0$ also for metal hydroxides. Add them to the theory-vs-experiment plot from the previous task. How accurate is the MACE prediction?
`````

## Solvated ions

For the solvated proton we'll use the expression for $\mu^0_\mathrm{H^+}$ you derived previously. For the aqueous ions, we'll use experimental energies (don't forget to add the logarithmic concentration term).

`````{admonition} Task
:class: tip

Complete the Pourbaix diagram using MACE energies for your assigned materials. How does it compare to the experimental diagram, and to the one from the Materials project?
`````

`````{admonition} Optional tasks
:class: tip

* How does the Pourbaix diagram of r2SCAN-MACE compare to PBE-MACE?

* How does PBE-MACE compare to other PBE models? You can find a ranking with links to code at https://matbench-discovery.materialsproject.org/.

* How does geometry optimization affect the Pourbaix diagram?
`````


## References

```{bibliography}
:style: unsrt
:filter: docname in docnames
```
