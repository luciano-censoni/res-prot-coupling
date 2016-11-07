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


def res_prot_coupling(lap_mat, length, k_r, tau):
  neg_lap_mat = (-1.0 * lap_mat)
  coupling = []
  T0 = matrix(zeros(length)) #matrix type for transposition
  bath_mat = zeros((length, length))
  for i in range(length):
    T0 *= 0.0       #resetting
    bath_mat *= 0.0 #resetting
    T0[0,i] = 300.0 #heated residue
    bath_mat[i][i] =  k_r #B matrix
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
  k_r = float(raw_input("Input a value for k_r: "))

try:
  tau = float(argv[3])
except:
  tau = float(raw_input("Input a value for tau: "))

print "Protein file:", fil
print "k_r =", k_r
print "tau =", tau

ident, _, ext = fil.rpartition(".")
del _
ident = ident.rpartition("/")[-1]

if ext == "pdb": #naive extension file type check
  adj_mat, resid_name = contacts(fil, hydrogens, cutoff)
else:
  print "Unknown file type; terminating."
  exit()
#reading pdb file into adjacency matrix

assert adj_mat.shape[0] == len(resid_name) #size
length = len(resid_name)

for i in range(length):
  adj_mat[i][i] = 0.0 #no loops

try:
  savetxt(ident+".matrix.dat", adj_mat, fmt="%5.3f")
  print "Adjacency matrix successfully written."
except:
  print "Problem writing adjacency matrix; terminating."
  exit()

if use_nx:
  graph = Graph()

  for i in range(length):
    graph.add_node( resid_name[i][0], residue = resid_name[i][1] + str(resid_name[i][0]) )
  #possibly different numbering stops creating graph directly from adjacency matrix

  pairs = []
  for i in range(length):
    for j in range(length):
      if adj_mat[i][j]: pairs.append( (resid_name[i][0],resid_name[j][0],{'weight':1.0/adj_mat[i][j]}) )
      #edge weight as a function of distance may be changed here
  graph.add_edges_from(pairs)

  try:
    write_gpickle(graph, ident+".graph.gpickle")
    print "Graph written successfully."
  except:
    print "Problem writing graph; terminating."
    exit()
#use_nx

#for debugging
#k_r = 10.0**3; tau = 14.5 #2VUJ
#k_r = 10.0**3; tau = 11.0 #2VUL
#k_r = 10.0**3; tau = 15.0 #1F5J
#k_r = 10.0**3; tau = 18.0 #1M4W
#k_r = 10.0**3; tau = 9.5  #1XNB

if use_nx:
  lap_mat = laplacian_matrix(graph).toarray() #original type is sparse matrix
else:
  lap_mat = (-1.0 * adj_mat)
  for i in range(length):
    lap_mat[i][i] = sum(adj_mat[i])

print "Calculating couplings..."
coupling = zscore( res_prot_coupling(lap_mat, length, k_r, tau) )

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

#checking for missing residues
for i in range(1,len(nod)):
  if nod[i]-nod[i-1] != 1:
    print "WARNING: possible missing residues in pdb file. Gap found between residues", nod[i-1], "and", nod[i]
    print  "Careful interpretation of all generated plots is recommended."

print "Data saved successfully."
print "Program end."
