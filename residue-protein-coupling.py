#AUTHOR: Luciano Censoni, 2016; luciano.censoni at gmail dot com
#Institute of Chemistry, University of Campinas
#USAGE: python residue-protein-coupling.py FILE.pdb k_r tau

from contacts import contacts
from sys import argv, exit
from numpy import array, savetxt, zeros, matrix
from scipy.linalg import expm
from scipy.stats import zscore
from matplotlib import pyplot as py

use_nx = False
#set to True to use the networkx package and write an additional graph 'gpickle' file
if use_nx:
  from networkx import Graph, write_gpickle, laplacian_matrix


def res_prot_coupling(lap_mat, length, tau, k_r):
  neg_lap_mat = -1.0 * lap_mat
  coupling = []
  for i in range(length):
    T0 = zeros(length)
    T0[i] = 1.0 #heated residue
    T0 = matrix(T0)
    bath_mat = (T0.T * T0) * k_r #B matrix
    T0 *= 300.0
    coupling.append( ( (expm( (neg_lap_mat - bath_mat) * tau) * (T0.T-300.0)) + 300.0 ).mean() )
  return coupling


hydrogens = False
#'hydrogens' flag to control whether to consider hydrogen atoms
#hydrogens == True --> consider hydrogens
#hydrogens == False -> ignore hydrogens

cutoff = 6.0

#main routine body begins
try:
  fil = argv[1]
except:
  fil = raw_input("Input protein .pdb file path: ")

try:
  k_r = float(argv[2])
except:
  k_r = raw_input("Input a value for k_r: ")

try:
  tau = float(argv[3])
except:
  tau = raw_input("Input a value for tau: ")

ident, _, ext = fil.rpartition(".")
del _
ident = ident.rpartition("/")[-1]

if ext == "pdb": #naive extension file type check
  adj_mat, resid_name = contacts(fil, hydrogens, cutoff)
else:
  print "Unknown file type; terminating."
  exit()
#reading pdb file into adjacency matrix

for i in range(adj_mat.shape[0]):
  adj_mat[i][i] = 0.0 #no loops

try:
  savetxt(ident+".matrix.dat", adj_mat, fmt="%5.3f")
  print "Adjacency matrix successfully written."
except:
  print "Problem writing adjacency matrix; terminating."
  exit()

assert adj_mat.shape[0] == len(resid_name) #size
length = len(resid_name)

if use_nx:
  graph = Graph()

  for i in range(length):
    graph.add_node( resid_name[i][0], residue = resid_name[i][1] + str(resid_name[i][0]) )
  #possibly different numbering stops creating graph directly from adjacency matrix

  pairs = []
  for i in range(length):
    for j in range(length):
      if adj_mat[i][j]: pairs.append( (resid_name[i][0],resid_name[j][0],{'weight':1.0/adj_mat[i][j]}) )

  graph.add_edges_from(pairs)

  try:
    write_gpickle(graph, ident+".graph.gpickle")
    print "Graph written successfully."
  except:
    print "Problem writing graph; terminating."
    exit()
#use_nx

#for debugging
#tau = 14.5; k_r = 10.0**3 #2VUJ
#tau = 11.0; k_r = 10.0**3 #2VUL
#tau = 15.0; k_r = 10.0**3 #1F5J
#tau = 18.0; k_r = 10.0**3 #1M4W
#tau = 9.5;  k_r = 10.0**3 #1XNB

if use_nx:
  lap_mat = laplacian_matrix(graph).toarray() #original type is sparse matrix
else:
  lap_mat = -1.0 * adj_mat
  for i in range(length):
    row_sum = sum(adj_mat[i])
    lap_mat[i][i] = row_sum

print "Calculating couplings..."
coupling = zscore( res_prot_coupling(lap_mat, length, tau, k_r) )

if use_nx:
  nod = graph.nodes()
else:
  nod = map(lambda x: x[0], resid_name)
#plotting, saving
py.plot(nod, coupling)
py.grid(True)
py.title("Residue-protein thermal coupling, "+ident)
py.ylabel("Z-Score")
py.xlabel("Residue")
py.savefig(ident+"-coupling.png")
py.clf()

try:
  savetxt(ident+"-coupling.dat", array([nod, coupling]).T, fmt="%5.3f")
except:
  print "Problem writing couplings; terminating."
  exit()

print "Data saved successfully."
print "Program end."
