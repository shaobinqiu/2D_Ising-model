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
    infs=[]#Sij sum_si
    all_inf=[]#include filename

    for dir1 in list1:
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

            a=str(numbers[0])+'    '+str(numbers[1])+'    '\
                +'\n'
            b=str(numbers[0])+'    '+str(numbers[1])+'    '\
                +'    '+path+'/'+file+'\n'
            infs.append(a)
            all_inf.append(b)
            i += 1
            print(i)



    f = open('BC.txt', 'w')
    for inf in infs:
        f.write(inf)
    f.close
    f1 = open('allinf_BC.txt', 'w')
    for inf1 in all_inf:
        f1.write(inf1)
    f1.close

list1 = traversalDir_FirstDir("./")
print(list1)

main(list1)
