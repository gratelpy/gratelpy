GraTeLPy
========

GratTeLPy (Graph Theoretic Analysis of Linear Stability) is a software tool for parameter-free, graph-theoretic linear stability analysis. 
Given a mechanism file that describes a chemical reaction network (CRN) of mass-action reactions, GraTelPy analyzes the provided mechanism and determines if it meets a necessary condition for multistability.

This necessary condition is the existence of critical fragments [1] which are necessary for a fold bifurcation thus permitting multistability.

Typical usage looks like this::

After following the installation instructions in `INSTALL`, invoke on the command prompt `python` and import GraTelPy:
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

GraTelPy comes with mechanism files that you can find in the `mechanisms` subdirectory.

There are only three scripts, provided with GraTelPy, that you need to analyze and enumerate all critical fragments in your mechanism:
- `gratelpy_fragment_server`: this script enumerates all fragments in your mechanism
- `gratelpy_subclient`: this script goes through the list of fragments provided by `gratelpy_fragment_server` fragment-by-fragment, enumerates all subgraphs per fragment, computes criticality for a given fragment, and returns all data to `gratelpy_fragment_server`
- `gratelpy_check_data`: once `gratelpy_subclient` finished the list of fragments, `gratelpy_fragment_server` saves all critical fragments (including their subgraphs) to `MECHANISM_NAME.vsg`. `gratelpy_check_data` is then used to process this data file, make certain that your data contains no errors, and output discovered critical fragments in a readable format.

Let us look at these three steps with a mechanism file provided by GraTelPy, `mechanisms/reversible_substrate_inhibition.txt`.

Open two terminal windows. In the first window, invoke the fragment server:
```bash
$ gratelpy_fragment_server mechanisms/reversible_substrate_inhibition.txt 4

rank of stoichiometry matrix = 3
numer of valid fragments: 17
[(('s3', 's2', 's1'), ('w5', 'w3', 'w5')), (('s3', 's2', 's1'), ('w5', 'w3', 'w3')), 
(('s3', 's2', 's1'), ('w5', 'w3', 'w1')), (('s3', 's2', 's1'), ('w4', 'w3', 'w5')), 
(('s3', 's2', 's1'), ('w4', 'w3', 'w3')), (('s3', 's2', 's1'), ('w4', 'w3', 'w1')), 
(('s3', 's2', 's4'), ('w5', 'w3', 'w6')), (('s3', 's2', 's4'), ('w4', 'w3', 'w6')), 
(('s3', 's1', 's4'), ('w5', 'w5', 'w6')), (('s3', 's1', 's4'), ('w5', 'w3', 'w6')), 
(('s3', 's1', 's4'), ('w5', 'w1', 'w6')), (('s3', 's1', 's4'), ('w4', 'w5', 'w6')), 
(('s3', 's1', 's4'), ('w4', 'w3', 'w6')), (('s3', 's1', 's4'), ('w4', 'w1', 'w6')), 
(('s2', 's1', 's4'), ('w3', 'w5', 'w6')), (('s2', 's1', 's4'), ('w3', 'w3', 'w6')), 
(('s2', 's1', 's4'), ('w3', 'w1', 'w6'))]
```
The two arguments `gratelpy_fragment_server` requires are the location of the mechanism file `mechanisms/reversible_substrate_inhibition.txt` and the number of chemical species found therein `4`.
`gratelpy_fragment_server` then outputs the rank of the stoichiometry matrix, the total number of fragments it discovered in the mechanism, and all fragments.

Note that currently GraTelPy uses its own naming scheme for chemical species and reactions where chemical species are denoted by `s...` and chemical reactions are denoted by `w...`.
However, GraTelPy also outputs a dictionary file that translates between the naming scheme you use in your mechanism file and its internal naming scheme -- more on this below.

At this point, `gratelpy_fragment_server` waits for a `gratelpy_subclient` instance to start processing the 17 fragments it discovered. So in your second terminal window run the command
```bash
$ gratelpy_subclient mechanisms/reversible_substrate_inhibition.txt 4
Got (('s3', 's2', 's1'), ('w5', 'w3', 'w5'))
Got (('s3', 's2', 's1'), ('w5', 'w3', 'w3'))
Got (('s3', 's2', 's1'), ('w5', 'w3', 'w1'))
Got (('s3', 's2', 's1'), ('w4', 'w3', 'w5'))
Got (('s3', 's2', 's1'), ('w4', 'w3', 'w3'))
Got (('s3', 's2', 's1'), ('w4', 'w3', 'w1'))
Got (('s3', 's2', 's4'), ('w5', 'w3', 'w6'))
Got (('s3', 's2', 's4'), ('w4', 'w3', 'w6'))
Got (('s3', 's1', 's4'), ('w5', 'w5', 'w6'))
Got (('s3', 's1', 's4'), ('w5', 'w3', 'w6'))
Got (('s3', 's1', 's4'), ('w5', 'w1', 'w6'))
Got (('s3', 's1', 's4'), ('w4', 'w5', 'w6'))
Got (('s3', 's1', 's4'), ('w4', 'w3', 'w6'))
Got (('s3', 's1', 's4'), ('w4', 'w1', 'w6'))
Got (('s2', 's1', 's4'), ('w3', 'w5', 'w6'))
Got (('s2', 's1', 's4'), ('w3', 'w3', 'w6'))
Got (('s2', 's1', 's4'), ('w3', 'w1', 'w6'))
All done
```
Note that `gratelpy_subclient` requires the same arguments as `gratelpy_fragment_server`.

Also note that if  `gratelpy_fragment_server` discovers a large number of fragments, you can accelerate completion of the list of fragments simply by starting multiple instances of `gratelpy_subclient` just as above.

Further note that you  `gratelpy_fragment_server` and one or multiple instances of `gratelpy_subclient` can be run on different computers.
If `gratelpy_fragment_server` runs on a computer different from `gratelpy_subclient`, all you need to do is invoke `gratelpy_subclient` with the name of the computer (`hostname`) that `gratelpy_fragment_server` runs on
```bash
gratelpy_subclient mechanisms/reversible_substrate_inhibition.txt 4 HOSTNAME
```
This allows you to easily deploy GraTelPy on a large cluster with hundreds of subclients.

Once all fragments have been analyzed, `gratelpy_fragment_server` writes two files of interest to the directory from which you invoked `gratelpy_fragment_server`: a `.vsg` and a `.dict` file.

To output all discovered critical fragments in a readable format we now use the `gratelpy_check_data` script.

For the above example, GraTelPy writes files `reversible_substrate_inhibition.vsg` and `reversible_substrate_inhibition.dict`.

```bash
$ gratelpy_check_data reversible_substrate_inhibition.vsg
loading your data from reversible_substrate_inhibition.vsg

Name of your data is reversible_substrate_inhibition
Order of fragments in your data: 3
1 critical fragments reported
1 unique critical fragments
check multiplicity of complexes and reactions
checked multiplicities of 3 subgraphs (of total 3 )
checking if all paths in every subgraph are witihin a cycle ...
3 subgraphs tested and 3 of these have all paths in cycles
1 fragments checked of total 1 fragments. k_s scores ok: 1
printing first 20 critical fragments

('s3', 's2', 's1')
('w4', 'w3', 'w5')

printing all critical fragments with subgraphs and K_S score to file reversible_substrate_inhibition_order_3.dat
```
The output of this invocation of `gratelpy_check_data` shows that there is exactly one critical fragment, that all your data is correct, and that your critical fragments together with all constituent subgraphs are written to `reversible_substrate_inhibition_order_3.dat` which you can now open in a simple text editor.

The above invocation outputs the discovered critical fragments with GraTelPy's internal naming scheme. To output the critical fragments with the naming scheme you use in your mechanism file, provide the dictionary file `.dict`

```bash
$ gratelpy_check_data reversible_substrate_inhibition.vsg reversible_substrate_inhibition.dict

loading your data from reversible_substrate_inhibition.vsg
Name of your data is reversible_substrate_inhibition
loading dictionary file from reversible_substrate_inhibition.dict
Order of fragments in your data: 3
1 critical fragments reported
1 unique critical fragments
check multiplicity of complexes and reactions
checked multiplicities of 3 subgraphs (of total 3 )
checking if all paths in every subgraph are witihin a cycle ...
3 subgraphs tested and 3 of these have all paths in cycles
1 fragments checked of total 1 fragments. k_s scores ok: 1
printing first 20 critical fragments

('A3', 'A2', 'A1')
('k4', 'k3', 'k5')

printing all critical fragments with subgraphs and K_S score to file reversible_substrate_inhibition_order_3.dat
```

[1]: Maya Mincheva and Marc R. Roussel (2007): Graph-theoretic methods for the analysis of chemical and biochemical networks. I. Multistability and oscillations in ordinary differential equation models. Journal of Mathematical Biology. Volume 55, Issue 1, pp 61-86.
dx.doi.org/10.1007/s00285-007-0099-1
