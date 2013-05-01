GraTeLPy
========

GratTeLPy (Graph Theoretic Analysis of Linear Stability) is a software tool for parameter-free, graph-theoretic linear stability analysis. 
Given a mechanism file that describes a chemical reaction network (CRN) of mass-action reactions, GraTelPy analyzes the provided mechanism and determines if it meets a necessary condition for multistability.

This necessary condition is the existence of critical fragments [1] which are necessary for a fold bifurcation thus permitting multistability.

Typical usage looks like this::

On the command prompt, invoke `python` and import GraTelPy:
```bash
Python 2.7.1 (r271:86832, Apr 12 2011, 16:15:16) 
[GCC 4.6.0 20110331 (Red Hat 4.6.0-2)] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> import gratelpy
```
Validate that you loaded GraTelPy correctly have the latest version installed
```bash
>>> gratelpy.get_version()
'0.1.0'
```

A typical mechanism file provided for GraTelPy looks as follows
```text
# Reversible Substrate Inhibition
# Example 4, Mincheva and Roussel 2007, [1]
# 4 complexes
[A1] -> ; k1
 -> [A1]; k2
[A1] + [A2] -> [A3] ; k3
[A3] -> [A2] ; k4
[A1] + [A3] -> [A4] ; k5
[A4] -> [A1] + [A3] ; k6
```
where `[NAME]` denotes a chemical species, the arrow `->` separates the reactant from the product side, and `;` separates the reaction from the corresponding rate constant `k`.
You are free to name your chemical species and rate constants as you please.
Also note that lines starting with `#` are comments that allow you to annotate your mechanism files.


[1]: Maya Mincheva and Marc R. Roussel (2007): Graph-theoretic methods for the analysis of chemical and biochemical networks. I. Multistability and oscillations in ordinary differential equation models. Journal of Mathematical Biology. Volume 55, Issue 1, pp 61-86.
dx.doi.org/10.1007/s00285-007-0099-1
