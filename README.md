GraTeLPy
========

## Requirements for GraTeLPy

GraTeLPy requires Python 2.6 or 2.7 and NetworkX 1.6 or above and runs on Windows, Mac OSX, and Linux.

We test GraTeLPy continually on a variety of platforms to ensure that it will work for you.
If you encounter any problems installing or using GraTeLPy please open an
[issue](https://github.com/gratelpy/gratelpy/issues) or [drop us a line](mailto:gratelpy@gmail.com).

## GraTeLPy Status

We test continually GraTeLPy on Linux (64-bit Ubuntu Linux 12.04) and Mac OS X 10.8.5.
A green status icon means that the most recent version of GraTeLPy ran correctly for that software configuration.

On **Linux** with both Python 2.6 and 2.7 with NetworkX 1.6, 1.7, 1.8, and the latest version of NetworkX:
[![Build Status](https://travis-ci.org/gratelpy/gratelpy.png?branch=master)](https://travis-ci.org/gratelpy/gratelpy)

On **Mac OS X** with Python 2.7 (2.6 testing coming soon hopefully) 
with NetworkX 1.6, 1.7, 1.8, and the latest version of NetworkX:
[![Build Status](https://travis-ci.org/gratelpy/gratelpy.png?branch=master.osx)](https://travis-ci.org/gratelpy/gratelpy)

Testing on **Windows** is done by hand and therefore less frequent. 
Our latest test of the current version of GraTeLPy with the current version
of the [Continuum Analytics Anaconda Python bundle](http://continuum.io/downloads) on Windows 8 worked correctly.

## Download GraTeLPy

We improve GraTeLPy continually and release all tested and working improvements directly on GitHub.
The latest working version of GraTeLPy can always be downloaded
[here](https://github.com/gratelpy/gratelpy/archive/master.zip).

From time to time we release numbered versions to mark specific improvements and developments.
A list of our numbered releases can be found
[here](https://github.com/gratelpy/gratelpy/releases).

The latest release of these numbered versions is also always avaiable on PyPI
[here](https://pypi.python.org/pypi/GraTeLPy).

## Who Uses GraTeLPy

* [Md. Ruhul Amin, Marc R Roussel (2013): Graph-Theoretic Analysis of a Model for the Coupling 
Between Photosynthesis and Photorespiration](http://www.nrcresearchpress.com/doi/abs/10.1139/cjc-2013-0315)

## Documentation Note

**Please Note**: At present the documentation we provide (in doc/) is
very rudimentary. Until we update our documentation, this README.md file
will act in place of proper documentation that we owe you.

If anything remains unclear or is confusing, please open an
[issue](https://github.com/gratelpy/gratelpy/issues) or
[drop us a line](mailto:gratelpy@gmail.com) and we will do our best
to help.

## Introduction

GraTeLPy (GRAph ThEoretic Analysis of Linear Stability in Python -- slightly contrived) is a software tool for parameter-free, graph-theoretic linear stability analysis. 
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

**Note**: this will also upgrade NetworkX if that library has been updated.

## Test Your GraTeLPy Setup

GraTeLPy provides two separate scripts that the user can invoke directly from their command line
to test their installation of GraTeLPy.
On Linux and Mac OSX the command line is the terminal while on Windows the command line is the console.

The script `gratelpy_test` runs a series of tests that ensure that the implementation of GraTeLPy
works as expected. For instance, we know the rank of the stoichiometric matrices of the mechanisms
provided with GraTeLPy (in `gratelpy/mechanisms`). `gratelpy_test` computes the rank of all of these
mechansims and tests if the computed values are equal to the expected values.

To invoke `gratelpy_test` simply invoke it from your command prompt:
<!--- Ignore 'bash': this is only for syntax highlighting on
GitHub.com, https://help.github.com/articles/github-flavored-markdown#syntax-highlighting --->
```bash
gratelpy_test
```

Depending on your versions of Python and NetworkX you should see output similar to the following.
The important line here is the last one which should show `OK`.
If you encounter any errors, please [report those](https://github.com/gratelpy/gratelpy/issues).
<!--- Ignore 'bash': this is only for syntax highlighting on
GitHub.com, https://help.github.com/articles/github-flavored-markdown#syntax-highlighting --->
```bash

GraTeLPy Copyright (C) 2013  Georg R Walther, Matthew Hartley.
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it under certain conditions.
For details visit https://github.com/gratelpy/gratelpy and read our LICENSE or contact 
the author at gratelpy@gmail.com.

Your Python version is 2.6.
Your NetworkX version is 1.6.

............
----------------------------------------------------------------------
Ran 12 tests in 0.018s
OK
```

A second script, `gratelpy_time`, tests how fast GraTeLPy runs on your system.
This script enumerates all critical fragments for a number of the mechanisms provided with GraTeLPy
on both one processor and multiple processors (corresponding to one or multiple GraTeLPy clients).

To invoke `gratelpy_time` call it from your command prompt:
<!--- Ignore 'bash': this is only for syntax highlighting on
GitHub.com, https://help.github.com/articles/github-flavored-markdown#syntax-highlighting --->
```bash
gratelpy_time
```

Depending on your version of Python, NetworkX and your operating system you should see output
similar to the following (note that timings will vary for your system):

```bash

GraTeLPy Copyright (C) 2013  Georg R Walther, Matthew Hartley.
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it under certain conditions.
For details visit https://github.com/gratelpy/gratelpy and read our LICENSE or contact 
the author at gratelpy@gmail.com.

Your Python version is 2.7
Your NetworkX version is 1.9.dev_20131206011125.

One client:

Analyzing /usr/lib/python2.7/site-packages/GraTeLPy-0.2.0-py2.7.egg/gratelpy/
mechanisms/reversible_substrate_inhibition.txt ...
run time 0.0888 seconds
run time 0.0514 seconds
run time 0.0478 seconds
median 0.0514 seconds
Analyzing /usr/lib/python2.7/site-packages/GraTeLPy-0.2.0-py2.7.egg/gratelpy/
mechanisms/glycolysis_mechanism.txt ...
run time 5.8973 seconds
run time 5.8835 seconds
run time 5.9257 seconds
median 5.8973 seconds
Analyzing /usr/lib/python2.7/site-packages/GraTeLPy-0.2.0-py2.7.egg/gratelpy/
mechanisms/cdc42_yeast.txt ...
run time 10.0858 seconds
run time 10.2031 seconds
run time 10.3910 seconds
median 10.2031 seconds

Two clients:

Analyzing /usr/lib/python2.7/site-packages/GraTeLPy-0.2.0-py2.7.egg/gratelpy/
mechanisms/cdc42_yeast.txt ...
run time 7.4444 seconds
run time 7.4505 seconds
run time 7.2978 seconds
median 7.4444 seconds
Analyzing /usr/lib/python2.7/site-packages/GraTeLPy-0.2.0-py2.7.egg/gratelpy/
mechanisms/single_layer_mapk_mechanism.txt ...
run time 10.9153 seconds
run time 11.4756 seconds
run time 10.8592 seconds
median 10.9153 seconds
```

The above was invoked on a Linux computer. Invoking the same script on Mac OSX or Windows will show different
path locations for the mechanism `.txt` files.

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

<!--- Ignore 'bash': this is only for syntax highlighting on
GitHub.com, https://help.github.com/articles/github-flavored-markdown#syntax-highlighting --->
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

The two arguments that `gratelpy_fragment_server` requires are the location of the mechanism file `mechanisms/reversible_substrate_inhibition.txt` and the number of chemical species found therein `4`.
`gratelpy_fragment_server` then outputs the rank of the stoichiometry matrix, the total number of fragments it discovered in the mechanism, and all fragments.

Note that currently GraTeLPy uses its own naming scheme for chemical species and reactions where chemical species are denoted by `s...` and chemical reactions are denoted by `w...`.
However, GraTeLPy also outputs a dictionary file that translates between the naming scheme you use in your mechanism file and its internal naming scheme -- more on this below.

At this point, `gratelpy_fragment_server` waits for a `gratelpy_subclient` instance to start processing the 17 fragments it discovered. So in your second terminal window run the command

<!--- Ignore 'bash': this is only for syntax highlighting on
GitHub.com, https://help.github.com/articles/github-flavored-markdown#syntax-highlighting --->
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

<!--- Ignore 'bash': this is only for syntax highlighting on
GitHub.com, https://help.github.com/articles/github-flavored-markdown#syntax-highlighting --->
```bash

$ gratelpy_subclient mechanisms/reversible_substrate_inhibition.txt 4 HOSTNAME
```

This allows you to easily deploy GraTeLPy on a large cluster with hundreds of subclients.

Once all fragments have been analyzed, `gratelpy_fragment_server` writes two files of interest to the directory from which you invoked `gratelpy_fragment_server`: a `.vsg` and a `.dict` file.

To output all discovered critical fragments in a readable format we now use the `gratelpy_check_data` script.

For the above example, GraTeLPy writes files `reversible_substrate_inhibition.vsg` and `reversible_substrate_inhibition.dict`.

<!--- Ignore 'bash': this is only for syntax highlighting on
GitHub.com, https://help.github.com/articles/github-flavored-markdown#syntax-highlighting --->
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

The above invocation outputs the discovered critical fragments with GraTeLPy's internal naming scheme. To output the critical fragments with the naming scheme you use in your mechanism file, provide the dictionary file `.dict`

<!--- Ignore 'bash': this is only for syntax highlighting on
GitHub.com, https://help.github.com/articles/github-flavored-markdown#syntax-highlighting --->
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
