from numpy import *
import numpy as np
import itertools

def from_string(content):
    # move empty line
    lines = [l for l in content.split('\n') if l.rstrip()]
    comment = lines[0]
    zoom = float(lines[1])
    lattice = np.around(np.array([[float(i) for i in line.split()]
                                                 for line in lines[2:5]]),
                                 decimals=6)
    if zoom < 0:
        # In vasp, a negative scale factor is treated as a volume. We need
        # to translate this to a proper lattice vector scaling.
        vol = abs(np.linalg.det(lattice))
        lattice *= (-zoom / vol) ** (1 / 3)
    else:
        lattice *= zoom
    #nsymbols = [Specie(s).Z for s in lines[5].split()]
    natoms = [int(i) for i in lines[6].split()]
    positions = np.around(np.array([[float(i) for i in line.split()[0:3]]
                                                 for line in lines[8:]]),
                          decimals=6)
    return lattice, positions, natoms

with open('./sc.vasp', "r") as f:
    content = f.read()

def get_s(atom_number, numbers):
    # check the spin direction
    if atom_number > sum(numbers)-numbers[-1]:
        s = -1
    elif atom_number > sum(numbers)-numbers[-1]-numbers[-2]:
        s = 1
    else:
        s = 0
    return s


latt, pos, numbers = from_string(content)
print(latt, pos, numbers)
print(type(latt),type(pos),sum(numbers))
list1=range(1,sum(numbers)+1)
# gennerate a repeat matrix
mt_rep=np.array([[a_rep, b_rep, c_rep] for a_rep in [-1,0,1]
                                            for b_rep in [-1,0,1]
                                                for c_rep in [-1,0,1]])
pos_rep = np.array([vect+atom for vect in mt_rep
                                            for atom in pos])
print(pos_rep)
atom_number1 = 0
S_ij = 0
for atom in pos:
    atom_number1 += 1
    # direct vector of this atom and all atoms
    radius = mat(pos_rep-tile(atom, (len(pos_rep),1)))*mat(latt)
    diss = list(np.linalg.norm(radiu) for radiu in radius)
    step = 0
    s_i = get_s(atom_number1, numbers)
    diss_min = min(i for i in diss if i>0)
    for dis in diss:
        step += 1
        if dis == diss_min:
            s_j = get_s(((step-1) % sum(numbers))+1, numbers)
            S_ij += s_i * s_j
            print(s_i, s_j, step)

print(S_ij)
f = open('s_ij_s_i.txt', 'r')
