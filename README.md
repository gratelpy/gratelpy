GraTeLPy
========

**Please Note: At present the documentation we provide (in doc/) is
very rudimentary. Until we update our documentation, this README.md file
will act in place of proper documentation that we owe you.

If anything remains unclear or is confusing, please open an
[issue](https://github.com/gratelpy/gratelpy/issues) or
[drop us a line](mailto:gratelpy@gmail.com) and we will do our best
to help.
**

## Introduction

GratTeLPy (Graph Theoretic Analysis of Linear Stability) is a software tool for parameter-free, graph-theoretic linear stability analysis. 
Given a mechanism file that describes a chemical reaction network (CRN) of mass-action reactions, GraTelPy analyzes the provided mechanism and determines if it meets a necessary condition for multistability.

This necessary condition is the existence of critical fragments [1] which are necessary for a fold bifurcation thus permitting multistability.

## Installation

To install GraTeLPy please make sure you have both Python (at least 2.7) 
and NetworkX (at least version 1.7) installed.

You can either install GraTeLPy from PyPI (the Python Package Index) 
directly if you have `pip` installed or
first download the source code and then run the setup on your computer.

### Installation from PyPI Using `pip`

If you have permission to write to root and if you wish to install
GraTeLPy system-wide you can invoke the following from your command line

    $ pip install gratelpy

This will most likely place the provided binaries in `/usr/local/bin` and
our Python source code in `/usr/lib/pythonX.Y/site-packages` (X.Y = 2.7 if
Python 2.7 is installed).

If you wish to install GraTeLPy locally to your home directory, invoke

    $ pip install gratelpy --user
       
This will likely place both the binaries and the source code in
`$HOME/.local/bin` and `$HOME/.local/lib/pythonX.Y/site-packages/`
respectively. Please ensure that your `PATH` environment variable contains
`$HOME/.local/bin` and that `PYTHONPATH` includes
`$HOME/.local/lib/pythonX.Y/site-packages/`.

### Installation from Source

To install directly from source, either clone into our git repository

    $ git clone https://github.com/gratelpy/gratelpy/archive/master.zip
    
or [download](https://github.com/gratelpy/gratelpy/archive/master.zip)
the current version of our code in a zip file.

Direct your terminal to the top-level directory of GraTeLPy and make certain
that you can locate the `setup.py` file.

For a system-wide installation (as above), invoke the following

    $ python setup.py install

and for a local installation (as above), append the `--user` option

    $ python setup.py install --user
    
### Where Are the GraTeLPy Binaries and Source Code?

Both processes outlined above will output where both the binaries and
the source code are copied to.
These locations may vary widely (as has been reported by users) and we are
afraid that at this point we do not fully understand what causes these
variations.

There are, however, ways to enforce that both setup methods (`pip` and
`python setup.py install`) copy all relevant files to a specific directory.

    $ pip install --install-option="--prefix=PATH" gratelpy
    $ python setup.py install --prefix=PATH

Where `PATH` is the *prefix* where you wish GraTeLPy to be installed.
If, for instance, you want GraTeLPy to be installed to `$HOME/MY_LOCAL`
instead of `$HOME/.local` then use either one of these methods:

    $ pip install --install-option="--prefix=$HOME/MY_LOCAL" gratelpy
    $ python setup.py install --prefix=$HOME/MY_LOCAL

The binaries and source code are then in `$HOME/MY_LOCAL/bin` and
`$HOME/MY_LOCAL/lib/pythonX.Y/site-packages` respectively.

Note that `python setup.py install` also gives you a *dry-run* option which
allows you to see where GraTeLPy will be installed to without carrying out
the actual copy operations:

    $ python setup.py install --dry-run
    $ python setup.py install --dry-run --user

Unfortunately, this option does not appear to be available for `pip`.

## Upgrading GraTeLPy

Use the mode of installation you chose from the above options and append the
`--upgrade` flag:

    $ pip install gratelpy --upgrade
    $ pip install gratelpy --user --upgrade
    $ python setup.py install --upgrade
    $ python setup.py install --user --upgrade

** Note: this will also upgrade NetworkX if that library has been updated **

## Test Your GraTeLPy Setup

Invoke on the command prompt `python` and import GraTeLPy:
```bash
<!--- Ignore 'bash': this is only for syntax highlighting on
GitHub.com, https://help.github.com/articles/github-flavored-markdown#syntax-highlighting --->

Python 2.7.1 (r271:86832, Apr 12 2011, 16:15:16) 
[GCC 4.6.0 20110331 (Red Hat 4.6.0-2)] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> import gratelpy
```
Validate that you loaded GraTeLPy correctly and that you have 
the latest version installed
```bash
<!--- Ignore 'bash': this is only for syntax highlighting on
GitHub.com, https://help.github.com/articles/github-flavored-markdown#syntax-highlighting --->

>>> gratelpy.get_version()
'0.1.1'
```

## Mechanism Files

A typical mechanism file for GraTeLPy looks as follows
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

GraTelPy comes with mechanism files that you can find in the `mechanisms/` subdirectory.

## Working with GraTeLPy Binaries / Scripts

There are only three scripts, provided with GraTelPy, that you need to analyze and enumerate all critical fragments in your mechanism:
- `gratelpy_fragment_server`: this script enumerates all fragments in your mechanism
- `gratelpy_subclient`: this script goes through the list of fragments provided by `gratelpy_fragment_server` fragment-by-fragment, enumerates all subgraphs per fragment, computes criticality for a given fragment, and returns all data to `gratelpy_fragment_server`
- `gratelpy_check_data`: once `gratelpy_subclient` finished the list of fragments, `gratelpy_fragment_server` saves all critical fragments (including their subgraphs) to `MECHANISM_NAME.vsg`. `gratelpy_check_data` is then used to process this data file, make certain that your data contains no errors, and output discovered critical fragments in a readable format.

Let us look at these three steps with a mechanism file provided by GraTeLPy, `mechanisms/reversible_substrate_inhibition.txt`.

Open two terminal windows. In the first window, invoke the fragment server:
```bash 
<!--- Ignore 'bash': this is only for syntax highlighting on
GitHub.com, https://help.github.com/articles/github-flavored-markdown#syntax-highlighting --->

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

The two arguments that `gratelpy_fragment_server` requires are the location of the mechanism file `mechanisms/reversible_substrate_inhibition.txt` and the number of chemical species found therein `4`.
`gratelpy_fragment_server` then outputs the rank of the stoichiometry matrix, the total number of fragments it discovered in the mechanism, and all fragments.

Note that currently GraTeLPy uses its own naming scheme for chemical species and reactions where chemical species are denoted by `s...` and chemical reactions are denoted by `w...`.
However, GraTeLPy also outputs a dictionary file that translates between the naming scheme you use in your mechanism file and its internal naming scheme -- more on this below.

At this point, `gratelpy_fragment_server` waits for a `gratelpy_subclient` instance to start processing the 17 fragments it discovered. So in your second terminal window run the command
```bash
<!--- Ignore 'bash': this is only for syntax highlighting on
GitHub.com, https://help.github.com/articles/github-flavored-markdown#syntax-highlighting --->

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
<!--- Ignore 'bash': this is only for syntax highlighting on
GitHub.com, https://help.github.com/articles/github-flavored-markdown#syntax-highlighting --->

$ gratelpy_subclient mechanisms/reversible_substrate_inhibition.txt 4 HOSTNAME
```
This allows you to easily deploy GraTeLPy on a large cluster with hundreds of subclients.

Once all fragments have been analyzed, `gratelpy_fragment_server` writes two files of interest to the directory from which you invoked `gratelpy_fragment_server`: a `.vsg` and a `.dict` file.

To output all discovered critical fragments in a readable format we now use the `gratelpy_check_data` script.

For the above example, GraTeLPy writes files `reversible_substrate_inhibition.vsg` and `reversible_substrate_inhibition.dict`.

```bash
<!--- Ignore 'bash': this is only for syntax highlighting on
GitHub.com, https://help.github.com/articles/github-flavored-markdown#syntax-highlighting --->

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

The above invocation outputs the discovered critical fragments with GraTeLPy's internal naming scheme. To output the critical fragments with the naming scheme you use in your mechanism file, provide the dictionary file `.dict`

```bash
<!--- Ignore 'bash': this is only for syntax highlighting on
GitHub.com, https://help.github.com/articles/github-flavored-markdown#syntax-highlighting --->

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
