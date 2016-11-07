# res-prot-coupling

Required packages: numpy, scipy, matplotlib
Optional networkx package may be found at https://networkx.readthedocs.io/en/stable/download.html

usage: python residue-protein-coupling.py FILE.pdb kr tau

(typical values are kr = 1000.0 and tau = 15.0)

The script will output the residue-protein relative thermal couplings in a ".dat" file, and use matplotlib to plot the couplings in a ".png" file. It will also write the adjacency matrix for the residue network in a ".dat" file.

If the networkx package is available, the script will write an additional "gpickle" graph file with the network specification to allow subsequent analysis. In order to use networkx please change the appropriate flag in line 12 of residue-protein-coupling.py.
