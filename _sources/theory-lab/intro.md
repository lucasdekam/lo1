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


```{note} 
:class: dropdown

Oxygen reduction can also take place through the binding of an entire oxygen molecule {cite:p}`norskov2004origin`, 

$$\mathrm{O_2 + * \to O_2^* }.$$

This is particularly relevant for porphyrins {cite:p}`lv2021controlling`.
```

The site $*$ acts as a catalyst: it binds the *intermediates* $\mathrm{O^*}$ and $\mathrm{HO^*}$ during the reaction. The reduction of oxygen to water causes a flow of electrons: an electric current. The field of science that studies catalytic reactions that drive or require an electric current is called *electrocatalysis*. 


## The binding energy in electrocatalysis

In the lab part, you will find or have found that different catalysts (different electrodes) have a different activity for the oxygen reduction reaction. A major milestone in electrocatalysis was relating the activity of a catalyst to the *binding energy* to an intermediate.

The binding energies of various oxygen species are related: if a catalyst binds oxygen strongly, it will also bind OH strongly {cite:p}`norskov2004origin, kulkarni2018understanding`.

For optimal catalysis, the catalyst should bind oxygen strongly enough to adsorb it on the surface and enable electron transfer, but not so strongly that the product (water) cannot leave the surface. 

Plotting the activity of various catalysts against their binding energy reveals a 'volcano plot': the activity is optimal if the binding energy is not too large, and not too small.


```{figure} ../images/volcano.png
---
height: 300
name: volcano
---
Volcano plot for the oxygen reduction reaction. 
```


**Task for students who model metals**: By writing down the Gibbs energy of the reaction for 

$$\mathrm{\frac12 O_2 + * \to O^*}$$

define the binding energy of an oxygen atom to a metal surface.


**Task for students who model porphyrins**: By writing down the Gibbs energy of the reaction for 

$$\mathrm{O_2 + * \to O_2^* }$$

define the binding energy of an oxygen molecule to a porphyrin catalyst.


## References

```{bibliography}
:style: unsrt
:filter: docname in docnames
```
