#AUTHOR: Luciano Censoni, 2016
#Institute of Chemistry - University of Campinas

from re import match
from numpy import array, dstack, logical_or, transpose
from scipy.spatial.distance import cdist
from os import getcwd
from json import dump
from itertools import compress


verbose = True
#verbose == True: all prints
#verbose == False: no prints


def get_atoms_pdb(pdb_file, hydrogens=False):
  #'hydrogens' flag to control whether to consider hydrogen atoms
  #hydrogens == True --> consider hydrogens
  #hydrogens == False -> ignore hydrogens

  try:
    pdb = open(pdb_file)
  except:
    raise IOError

  residues = [] #atom-residue index reference
  n_residues = 0 #number of residues
  current_res = None #current residue
  n_atoms = 0 #number of atoms

  pdb_re = '(.{6})(.{5})(.{5})(.{1})(.{4})(.{1})(.{4})(.{1})(.{11})(.{8})(.{8})(.{6})(.{6})(.{12})(.{2})'
  #regular expression, fields in an ATOM line in a pdb file. Must sum to 80 chars.

  mask = [] #boolean; False for atoms we want to ignore

  x = []
  y = []
  z = []
  #atomic positions

  resid_name = []
  #residue names and indexes

  for line in pdb:
    if line[0:3] == 'TER': break #chain end, exit loop
    if line[0:4] != 'ATOM' and line[0:6] != 'HETATM': continue #protein atoms not reached yet, jump to next line

    line = line[:-1] #removing \n from end of line
    while len(line) < 80: line = line + ' ' #charge field usually empty, fill with whitespace if needed

    l = match(pdb_re, line).groups() #separating fields

    if l[4] == 'TIP3': break #for pdbs with waters before TER
    if l[3] == 'B': continue #ignoring alternative conformations

    #having made it this far, this is a protein atom
    n_atoms += 1 #for consistency checking
    mask.append(True)

    if not hydrogens and l[13][-1] == 'H': mask[-1] = False #(not hydrogens == True) when hydrogens must be ignored
    #different conditions may be applied here

    #generating atom-residue reference:
    if current_res != int(l[6]): #first atom of first residue OR new residue has begun
      current_res = int(l[6])
      n_residues += 1 #also for consistency checking
      resid_name.append( [ int(l[6]), l[4].strip() ] ) #pdb residue number, name (first residue possibly not 1)

    residues.append(current_res) #atom-residue index reference

    #read coordinates
    x.append(float(l[8]))
    y.append(float(l[9]))
    z.append(float(l[10]))
  #end of atom loop

  pdb.close()

  if verbose: #for debugging
    print 'Number of atoms:', n_atoms
    print 'Number of residues:', n_residues

  #sanity checks
  assert len(x) == len(y)
  assert len(y) == len(z)
  assert len(z) == n_atoms
  assert n_atoms == len(mask)
  assert len(mask) == len(residues)

  assert len(resid_name) == n_residues

  x = list(compress(x, mask))
  y = list(compress(y, mask))
  z = list(compress(z, mask))
  residues = list(compress(residues, mask))

  n_atoms = len(x) #updating after masking

  #obtain list with the index of the first atom of each residue from residues, after masking
  indexes = [0]
  current = residues[0]
  for i in range(len(residues)):
    if current != residues[i]:
      indexes.append(i)
      current = residues[i]
  del residues

  return (x, y, z, n_atoms, n_residues, resid_name, indexes)
#get_atoms_pdb


def resid_adj_mat(dist_mat, cutoff, indexes):
  #indexes contains the index of the first atom of each residue
  if verbose: print 'Extracting contacts from distance matrix...'

  atom_adj_mat = (dist_mat <= cutoff)

  f = lambda x: logical_or.reduceat(x, indexes, axis=0)
  g = lambda x: logical_or.reduceat(x, indexes, axis=1)
  resid_adj_mat = f(g(atom_adj_mat)).astype(float)

  return resid_adj_mat
#resid_adj_mat


def contacts(pdb_file, hydrogens, cutoff):
  #'hydrogens' flag to control whether to consider hydrogen atoms
  #hydrogens == True --> consider hydrogens
  #hydrogens == False -> ignore hydrogens

  #read pdb
  x, y, z, n_atoms, n_residues, resid_name, indexes = get_atoms_pdb(pdb_file, hydrogens)

  #stack dimensions for cdist call
  pos = dstack((x, y, z))[0]

  dist_mat = cdist(pos, pos) #distances between all pairs of atoms

  return resid_adj_mat(dist_mat, cutoff, indexes), resid_name
#contacts
