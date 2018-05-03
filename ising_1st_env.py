from numpy import *
import numpy as np
import itertools
import os
from scipy.spatial.distance import pdist

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

#定义一个函数，path为你的路径
def traversalDir_FirstDir(path):
#定义一个列表，用来存储结果
    list1 = []
#判断路径是否存在
    if (os.path.exists(path)):
    #获取该目录下的所有文件或文件夹目录
        files = os.listdir(path)
        for file1 in files:
            #得到该文件下所有目录的路径
            m = os.path.join(path,file1)
            #判断该路径下是否是文件夹
            if (os.path.isdir(m)):
                h = os.path.split(m)
                list1.append(h[1])
    return list1

def main(list1):
    sub_str=np.array([[1, -1, -1, -1, -1, -1],[1, 1, -1, -1, -1, -1],[1, -1, 1, -1, -1, -1],
                       [1, -1, -1, 1, -1, -1],[1, 1, 1, -1, -1, -1],[1, -1, 1, 1, -1, -1],
                       [1, -1, 1, -1, 1, -1],[1, 1, 1, 1, -1, -1],[1, 1, 1, -1, 1, -1],
                       [1, 1, -1, 1, 1, -1],[1, 1, 1, 1, 1, -1],[1, 1, 1, 1, 1, 1],[-1, -1, 1, 1, -1, 1]])
    infs=[]#Sij sum_si
    all_inf=[]#include filename
    diss_min1st = 1.60
    diss_min2st = 2.77
    diss_min3st = 3.20
    diss_min4st = 4.233
    diss_min5st = 4.80
    diss_min6st =  5.54
    for dir1 in list1:
        R=3#扩胞半径
        path1 = "./" #文件夹目录
        path2 = dir1
        path =  path1 +path2
        files= os.listdir(path) #得到文件夹下的所有文件名称

        i=1
        for file in files: #遍历文件夹

            with open(path+'/'+file, "r") as f:
                content = f.read()
            latt, pos, numbers = from_string(content)
            if len(numbers)==1:
                numbers.append(0)
            number_atom= sum(numbers)
            list1=range(1,number_atom+1)
            # gennerate a repeat matrix
            costheta = np.dot(np.array(latt)[0], np.array(latt)[1])\
                         /(math.sqrt((np.array(latt)[0]**2).sum())*math.sqrt((np.array(latt)[1]**2).sum()))
            L = R/math.sqrt(1-costheta**2)
            n1 = round(L/math.sqrt((np.array(latt)[0]**2).sum()))+1
            n2 = round(L/math.sqrt((np.array(latt)[1]**2).sum()))+1
            mt_rep=np.array([[a_rep, b_rep, c_rep] for a_rep in range(-n1, n1+1)
                                                         for b_rep in range(-n2, n2+1)
                                                             for c_rep in [0]])

            pos_rep = np.array([vect+atom for vect in mt_rep
                                                        for atom in pos])
            atom_number1 = 0
            S_ij = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            for atom in pos:
                S = [0, 0, 0, 0, 0, 0]
                atom_number1 += 1
                # direct vector of this atom and all atoms
                radius = np.round(mat(pos_rep-tile(atom, (len(pos_rep),1)))*mat(latt),4)
                diss = list(math.sqrt((np.array(radiu)**2).sum()) for radiu in radius)

                #diss = list(np.linalg.norm(radiu) for radiu in radius)
                step = 0
                s_i = 2*bool(atom_number1 <=numbers[0])-1
                number_check = 0
                for dis in diss:
                    step += 1
                    if abs(dis-diss_min1st)<0.1:
                        s_j = s_i*(2*bool(((step-1) % number_atom)+1 <=numbers[0])-1)
                        if number_check==0:
                            S[0] = s_j
                            vector1 = step-1
                        else:
                            print(radius[step-1][:],radius[vector1][:])
                            cos = (radius[step-1][0]*radius[vector1][0]+radius[step-1][1]*radius[vector1][1])/dis**2
                            sin = np.sign(radius[vector1][0]*radius[step-1][1]-radius[vector1][1]*radius[step-1][0])*sqrt(1-(abs(cos)-0.001)**2)
                            #area = (bool(cos>=0.5)+2*bool(-0.5<=cos<0.5)+3*bool(cos<-0.5))*np.sign(sin)-bool(sin>0)
                            area = (bool(abs(cos-1)<0.05 and abs(sin-0)<0.05)+2*bool(abs(cos-0.5)<0.05 and abs(sin-0.866)<0.05)
                              +3*bool(abs(cos+0.5)<0.05 and abs(sin-0.866)<0.05)+4*bool(abs(cos+1)<0.05 and abs(sin-0)<0.05)
                              +5*bool(abs(cos+0.5)<0.05 and abs(sin+0.866)<0.05)+6*bool(abs(cos-0.5)<0.05 and abs(sin+0.866)<0.05))
                            S[area-1]=s_j
                            print(area,cos,sin)
                        number_check += 1

                if number_check != 6:

                    print(path+'/'+file+'\n')
                    baochuo

                for r in range(0,13):
                    for s in range(0,6):
                        if abs(sub_str[r][1]-S[1-s])+abs(sub_str[r][2]-S[2-s])+abs(sub_str[r][3]-S[3-s])+abs(sub_str[r][4]-S[4-s])+abs(sub_str[r][5]-S[5-s])+abs(sub_str[r][0]-S[0-s])==0:
                            if r != 12:
                                S_ij[r] += 1
                                print("C")
                            else:
                                S_ij[5] += 1
                                print("C")
                            break
                print(S)
            print(S_ij)
            a=str(S_ij[0]/number_atom)+'    '+str(S_ij[1]/number_atom)+'    '\
              +str(S_ij[2]/number_atom)+'    '+str(S_ij[3]/number_atom)+'    '\
              +str(S_ij[4]/number_atom)+'    '+str(S_ij[5]/number_atom)+'    '\
              + str(S_ij[6]/number_atom)+'    '+str(S_ij[7]/number_atom)+'    '\
                +str(S_ij[8]/number_atom)+'    '+str(S_ij[9]/number_atom)+'    '\
                +str(S_ij[10]/number_atom)+'    '+str(S_ij[11]/number_atom)+'    '\
                +'\n'
            b=str(S_ij[0]/number_atom)+'    '+str(S_ij[1]/number_atom)+'    '\
              +str(S_ij[2]/number_atom)+'    '+str(S_ij[3]/number_atom)+'    '\
              +str(S_ij[4]/number_atom)+'    '+str(S_ij[5]/number_atom)+'    '\
              + str(S_ij[6]/number_atom)+'    '+str(S_ij[7]/number_atom)+'    '\
                +str(S_ij[8]/number_atom)+'    '+str(S_ij[9]/number_atom)+'    '\
                +str(S_ij[10]/number_atom)+'    '+str(S_ij[11]/number_atom)+'    '\
                +'    '+path+'/'+file+'\n'
            infs.append(a)
            all_inf.append(b)
            i += 1
            print(i)



    f = open('infs_ev.txt', 'w')
    for inf in infs:
        f.write(inf)
    f.close
    f1 = open('all_inf_ev.txt', 'w')
    for inf1 in all_inf:
        f1.write(inf1)
    f1.close

list1 = traversalDir_FirstDir("./")
print(list1)

main(list1)
