from numpy import *
import numpy as np
import itertools
import os

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
    infs=[]#Sij sum_si
    all_inf=[]#include filename
    for dir1 in list1:
        R=6#扩胞半径
        path1 = "./" #文件夹目录
        path2 = dir1
        path =  path1 +path2
        files= os.listdir(path) #得到文件夹下的所有文件名称

        process=1#在算第几个结构
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
                                                             for c_rep in [0]])#扩胞矩阵
            print(mt_rep)


            s=[]
            for i in range(1,sum(numbers)*len(mt_rep)+1):
                s.append(1*bool(i>sum(numbers)-numbers[-1]-numbers[-2] and i<=sum(numbers)-numbers[-1])\
                    -1*bool(i>sum(numbers)-numbers[-1]))
            pos_rep = np.array([vect+atom for vect in mt_rep
                                                        for atom in pos])
            s_rep = tile(s, (1,len(pos_rep)))
            atom_number1 = 0
            print(len(s))
            S_ij = [0, 0, 0, 0, 0, 0]
            for atom in pos:
                atom_number1 += 1
                # direct vector of this atom and all atoms
                radius = mat(pos_rep-tile(atom, (len(pos_rep),1)))*mat(latt)
                diss = list(math.sqrt((np.array(radiu)**2).sum()) for radiu in radius)

                #diss = list(np.linalg.norm(radiu) for radiu in radius)
                step = 0
                s_i = s[atom_number1-1]
                diss_min1st = min(i for i in diss if i>0)
                diss_min2st = min(i for i in diss if i>diss_min1st+0.01)
                diss_min3st = min(i for i in diss if i>diss_min2st+0.01)
                diss_min4st = min(i for i in diss if i>diss_min3st+0.01)
                diss_min5st = min(i for i in diss if i>diss_min4st+0.01)
                diss_min6st = min(i for i in diss if i>diss_min5st+0.01)
                number_check = 0
                for dis in diss:
                    step += 1
                    # only consider nearst atoms
                    ord_dis=1*bool(abs(dis-diss_min1st)<0.001)+2*bool(abs(dis-diss_min2st)<0.001)\
                            +3*bool(abs(dis-diss_min3st)<0.001)+4*bool(abs(dis-diss_min4st)<0.001)\
                            +5*bool(abs(dis-diss_min5st)<0.001)+6*bool(abs(dis-diss_min6st)<0.001)
                    if ord_dis>0:
                        number_check += 1
                        s_j = s[step-1]
                        S_ij[ord_dis-1] += s_i * s_j
                print(number_check)
                if number_check % 6 != 0:
                    asjh
            a=str(S_ij[0]/number_atom)+'    '+str(S_ij[1]/number_atom)+'    '\
              +str(S_ij[2]/number_atom)+'    '+str(S_ij[3]/number_atom)+'    '\
              +str(S_ij[4]/number_atom)+'    '+str(S_ij[5]/number_atom)+'    '\
              +str((numbers[-2]-numbers[-1])/number_atom)+'    '+str(diss_min1st)+'    '\
              +str(diss_min2st)+'    '+str(diss_min3st)+'    '+str(diss_min4st)+\
              '    '+str(diss_min5st)+'    '+str(diss_min6st)+'\n'
            b=str(S_ij[0]/number_atom)+'    '+str(S_ij[1]/number_atom)+'    '\
              +str(S_ij[2]/number_atom)+'    '+str(S_ij[3]/number_atom)+'    '\
              +str(S_ij[4]/number_atom)+'    '+str(S_ij[5]/number_atom)+'    '\
              +str((numbers[-2]-numbers[-1])/number_atom)+'    '+str(diss_min1st)+'    '\
              +str(diss_min2st)+'    '+str(diss_min3st)+'    '+str(diss_min4st)+\
              '    '+str(diss_min5st)+'    '+str(diss_min6st)+'    '+path+'/'+file+'\n'
            infs.append(a)
            all_inf.append(b)
            process += 1
            print(process)



    f = open('infs.txt', 'w')
    for inf in infs:
        f.write(inf)
    f.close
    f1 = open('all_inf.txt', 'w')
    for inf1 in all_inf:
        f1.write(inf1)
    f1.close

list1 = traversalDir_FirstDir("./")
print(list1)

main(["tex"])
