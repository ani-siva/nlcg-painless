import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import math
import csv

def contour(X,Y,Z,X1,X2,X3,X4):
    fig,ax = plt.subplots(1,1)
    ax.contourf(Y,X,Z)
    ax.set_title("Contour Plot")
    ax.set_xlabel("X1 Variable")
    ax.set_ylabel("X2 Variable")
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.plot(X1,X2,'bo-',label='Gram Schmidt CG')
    plt.plot(X3,X4,'r.-',label='CG Method')
    plt.legend(loc='upper left')
    plt.show()

x1=[]
x2=[]
x3=[]
x4=[]
x5=[]
x6=[]

x=[]
y=[]
z=[]

with open("qf-output.txt") as textfile:
    reader=csv.reader(textfile,delimiter=' ')
    for row in reader:
        x.append(float(row[1]))
        y.append(float(row[2]))
        z.append(float(row[0]))

with open("gso-coords.txt") as textfile:
    reader=csv.reader(textfile,delimiter=' ')
    for row in reader:
        x1.append(float(row[0]))
        x2.append(float(row[1]))

with open("cg-method-coords.txt") as textfile:
    reader=csv.reader(textfile,delimiter=' ')
    for row in reader:
        x3.append(float(row[0]))
        x4.append(float(row[1]))

xi=np.unique(x)
yi=np.unique(y)
zi=np.array(z)

X,Y=np.meshgrid(xi,yi)
Z=np.reshape(zi,(np.size(X[0]),np.size(Y[0])))

X1=np.array(x1)
X2=np.array(x2)
X3=np.array(x3)
X4=np.array(x4)
X5=np.array(x5)
X6=np.array(x6)

contour(X,Y,Z,X1,X2,X3,X4)
