# Introduction

When you do an electrochemical experiment, the composition of your electrode can change depending on the conditions. For example, a metal electrode can oxidize, forming a metal oxide. Under different conditions, the metal can dissolve. The two main parameters affecting stability in a water-based electrolyte are **electrode potential** and **pH**. A Pourbaix diagram visualizes which material phase is stable under which conditions.

```{figure} ../../images/iron-pourbaix.png
---
height: 300
name: iron-pourbaix
---
Pourbaix diagram of iron ([Wikipedia](https://en.wikipedia.org/wiki/Pourbaix_diagram)).
```

## Water stability

The dashed blue lines in the Pourbaix diagram above indicate the water stability window. If you go below the lower line, you produce hydrogen according to the reaction

$$
    \mathrm{2H^+ + 2e^- \to H_2(g)}.
$$ (eq:her)

This reaction is called the **hydrogen evolution reaction** (HER). Above the upper line, the **oxygen evolution reaction** (OER)

$$
    \mathrm{2H_2O \to O_2(g) + 4H^+ + 4e^-}
$$ (eq:oer)

occurs. Typically, an experiment is done within the water stability window to avoid the formation of gas bubbles.

`````{admonition} Questions
:class: tip

* What iron oxide phases can you encounter in a typical electrochemical experiment?
* What do the HER and OER look like in alkaline media?

`````

## Reference points and SHE

Throughout this lab course, the concept of a 'reference point' will return several times. We cannot measure an absolute value for electric potentials and energies; we always measure the difference with respect to some reference point.

In electrochemistry, the most common reference point for the electrode potential is the **standard hydrogen electrode** (SHE). It's defined as a platinum electrode where the hydrogen evolution reaction {eq}`eq:her` occurs under standard conditions (1 bar $\mathrm{H_2}$-pressure, 1 M $\mathrm{H^+}$-concentration, i.e. pH 0). When using the SHE as a reference, the equilibrium potential of the SHE is defined as zero.

`````{admonition} Questions
:class: tip

* What point in the Pourbaix diagram corresponds to the SHE?
`````