from numpy import *
import numpy as np
import itertools
import os

import sys, time

class ProgressBar:
    # show rate of progress
    def __init__(self, count = 0, total = 0, width = 50):
        self.count = count
        self.total = total
        self.width = width
    def move(self):
        self.count += 1
    def log(self, s):
        sys.stdout.write(' ' * (self.width + 9) + '\r')
        sys.stdout.flush()
        progress = self.width * self.count / self.total
        sys.stdout.write('{0:3}/{1:3}: '.format(self.count, self.total))
        sys.stdout.write('#' * int(progress) + '-' * int(self.width - progress) + '\r')
        if progress == self.width:
            sys.stdout.write('\n')
        sys.stdout.flush()

bar = ProgressBar(total = 10)

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

def get_s(atom_number, numbers):
    # check the spin direction
    if atom_number > sum(numbers)-numbers[-1]:
        s = -1
    elif atom_number > sum(numbers)-numbers[-1]-numbers[-2]:
        s = 1
    else:
        s = 0
    return s
path1 = "/home/qiusb/Documents/python_work/2D_Ising_model/" #文件夹目录
path2 = "a"
path = path1 + path2
files= os.listdir(path) #得到文件夹下的所有文件名称
infs=[]#Sij sum_si
bar = ProgressBar(total = len(files))
i=1
for file in files: #遍历文件夹

    with open(path+'/'+file, "r") as f:
        content = f.read()
    latt, pos, numbers = from_string(content)
    list1=range(1,sum(numbers)+1)
    # gennerate a repeat matrix
    mt_rep=np.array([[a_rep, b_rep, c_rep] for a_rep in [-1,0,1]
                                                for b_rep in [-1,0,1]
                                                    for c_rep in [-1,0,1]])
    pos_rep = np.array([vect+atom for vect in mt_rep
                                                for atom in pos])
    atom_number1 = 0
    S_ij = 0
    for atom in pos:
        atom_number1 += 1
        # direct vector of this atom and all atoms
        radius = mat(pos_rep-tile(atom, (len(pos_rep),1)))*mat(latt)
        diss = list(np.linalg.norm(radiu) for radiu in radius)
        step = 0
        s_i = get_s(atom_number1, numbers)
        diss_min1st = min(i for i in diss if i>0)
        diss_min2st = min(i for i in diss if i>diss_min1st+0.1)
        diss_min3st = min(i for i in diss if i>diss_min2st+0.1)
        diss_min4st = min(i for i in diss if i>diss_min3st+0.1)
        diss_min5st = min(i for i in diss if i>diss_min4st+0.1)
        diss_min6st = min(i for i in diss if i>diss_min5st+0.1)
        print(diss,diss_min1st, diss_min2st, diss_min3st, diss_min4st, diss_min5st, diss_min6st)
        for dis in diss:
            step += 1
            # only consider nearst atoms
            if abs(dis-diss_min1st)<0.1:
                s_j = get_s(((step-1) % sum(numbers))+1, numbers)
                S_ij += s_i * s_j
    a=str(S_ij/sum(numbers))+'    '+str((numbers[-2]-numbers[-1])/sum(numbers))+'\n'
    infs.append(a)
    bar.move()
    bar.log('We have arrived at: ' + str(i+1))
    time.sleep(1)
    i += 1



f = open(path2+'_s_ij_s_i.txt', 'w')
for inf in infs:
    f.write(inf)
f.close
