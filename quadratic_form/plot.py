import matplotlib.pyplot as plt
import mpl_toolkits
from mpl_toolkits.mplot3d import Axes3D
import math
import numpy as np
import csv

def surfaceplot(X,Y,Z):
    my_cmap=plt.get_cmap('hot')
    ax=plt.axes(projection='3d')
    fig=plt.figure()
    surf=ax.plot_surface(X,Y,Z,cmap="plasma")
    cset=ax.contourf(X,Y,Z,zdir='z',offset=np.min(Z),cmap=my_cmap)
    cset=ax.contourf(X,Y,Z,zdir='x',offset=-5,cmap=my_cmap)
    cset=ax.contourf(X,Y,Z,zdir='y',offset=5,cmap=my_cmap)
    fig.colorbar(surf,ax=ax,shrink=0.3,aspect=5)
    ax.set_title('Surface and Contour Plot of the Quadratic Form Function')
    ax.set_xlabel('X1-variable')
    ax.set_xlim(-25,25)
    ax.set_ylabel('X2-variable')
    ax.set_ylim(-25,25)
    ax.set_zlabel('F(x) value')
    ax.set_zlim(np.min(Z),np.max(Z))
    plt.show()

def contour_plot(X,Y,Z):
    fig,ax=plt.subplots(1,1)
    cont=ax.contour(X,Y,Z)
    ax.set_title('Contour Plot of Quadratic Function')
    ax.set_xlabel('X1 variable')
    ax.set_ylabel('X2 variable')
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.show()

def surface():

    x=[]
    y=[]
    z=[]

    with open("qf-output.txt") as textfile:
        reader=csv.reader(textfile,delimiter=' ')
        for row in reader:
            x.append(float(row[1]))
            y.append(float(row[2]))
            z.append(float(row[0]))

    xi=np.unique(x)
    yi=np.unique(y)
    zi=np.array(z)

    X, Y=np.meshgrid(xi,yi)
    Z=np.reshape(zi,(np.size(X[0]),np.size(Y[0])))
    contour_plot(X,Y,Z)
    surfaceplot(X,Y,Z)

def vector_plot():
    x=[]
    y=[]
    zx=[]
    zy=[]

    with open("qf-vect-output.txt") as datafile:
        reader=csv.reader(datafile,delimiter=' ')
        for row in reader:
            x.append(float(row[2]))
            y.append(float(row[3]))
            zx.append(float(row[0]))
            zy.append(float(row[1]))
    xi=np.unique(x)
    yi=np.unique(y)
    zxi=np.array(zx)
    zyi=np.array(zy)
    X,Y=np.meshgrid(xi,yi)
    ZX=np.reshape(zxi,(np.size(X[0]),np.size(Y[0])))
    ZY=np.reshape(zyi,(np.size(X[0]),np.size(Y[0])))
    fig=plt.figure()
    ax=fig.add_subplot(1,1,1)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.set_xlabel('X1-Variable',fontweight='bold',loc="right")
    ax.set_ylabel('X2-Variable',fontweight='bold',loc="top")
    ax.set_title('Gradient Vector Field')
    plt.quiver(X,Y,ZX,ZY)
    plt.show()


surface()
