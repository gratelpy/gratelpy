.. hihglight::rst

Documentation
=============
Modules
-------
See link at top right for index of all modules provided by GraTelPy.

Mechanism Files
---------------
Mechanism files are written in plain text as standard chemical reactions -- one reaction per row.
Complexes (chemical substances) are encased in brackets, '[]', and each reaction is terminated with a semicolon, ';', followed by the rate constant.
Lines that start with the hash key, '#', are ignored when reading mechanism files hence can be used to add comments.

Example
^^^^^^^
|  # Reversible Substrate Inhibition
|  # Example 4, Mincheva and Roussel 2007
|  # 4 complexes
|  [A1] -> ; k1
|  -> [A1]; k2
|  [A1] + [A2] -> [A3] ; k3
|  [A3] -> [A2] ; k4
|  [A1] + [A3] -> [A4] ; k5
|  [A4] -> [A1] + [A3] ; k6
