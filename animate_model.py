import numpy as np
import glob 
from PIL import Image
import matplotlib.pyplot as plt
from scipy.stats.mstats import mquantiles
from mpl_toolkits.mplot3d import Axes3D
import sys
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import matplotlib.animation as animation

args = sys.argv
file_name = args[1]

LJ_CONVERTER = pow(2,1/6)

C_SIGMA = 3.851/LJ_CONVERTER
C_EPSILON = 0.105

O_SIGMA = 3.500/LJ_CONVERTER
O_EPSILON = 0.060

AU_SIGMA = 2.951
AU_EPSILON = 5.29

BR_SIGMA = 4.189/LJ_CONVERTER
BR_EPSILON = 0.251

S_SIGMA = 4.035/LJ_CONVERTER
S_EPSILON = 0.274

H_SIGMA = 2.886/LJ_CONVERTER
H_EPSILON = 0.044

ATOMS = {
    1 : [H_SIGMA, H_EPSILON, "white"],
    8 : [O_SIGMA, O_EPSILON, "blue"],
    12 : [C_SIGMA, C_EPSILON, "black"],
    16: [S_SIGMA, S_EPSILON, "yellow"],
    35 : [BR_SIGMA, BR_EPSILON, "darkred"],
    79 : [AU_SIGMA, AU_EPSILON, "gold"]
}

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.view_init(elev=0, azim=0)
ns = int(input())
np = int(input())
r = float(input())
xmax = 40
xmin = -40
ymax = 40
ymin = -40
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_zlim(-40, 40)
ax.set_box_aspect((80,80,80))
ppi=72 # points per inche 
ax_length=ax.bbox.get_points()[1][0]-ax.bbox.get_points()[0][0]
ax_point = ax_length*ppi/fig.dpi
xsize=xmax-xmin
fact=ax_point/xsize
X = []
Y = []
Z = []
N = []
for i in range(np):
    x, y, z, n = map(float, input().split())
    if z<=20 :
        n = int(n)
        X.append(x)
        Y.append(y)
        Z.append(z)
        N.append(n)
np = len(X)      
for i in range(ns):
    x, y, z, n = map(float, input().split())
    n = int(n)
    X.append(x)
    Y.append(y)
    Z.append(z)
    N.append(n)
    
def init():
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_zlim(-40, 40)
    ax.set_box_aspect((80,80,80))
    ppi=72 # points per inche 
    ax_length=ax.bbox.get_points()[1][0]-ax.bbox.get_points()[0][0]
    ax_point = ax_length*ppi/fig.dpi
    for i in range(np):
        x = X[i]
        y = Y[i]
        z = Z[i]
        n = N[i]
        ax.scatter(x, y, z, color=ATOMS[n][2],edgecolor='black',s=(fact*ATOMS[n][0])**2)
    for i in range(ns):
        x = X[i+np]
        y = Y[i+np]
        z = Z[i+np]
        n = N[i+np]
        ax.scatter(x, y, z, color=ATOMS[n][2],edgecolor='black',s=(fact*ATOMS[n][0])**2)
    return fig

def animate(i):
    ax.view_init(elev=-5,azim=3.6*i)
    return fig

ani = animation.FuncAnimation(fig, animate, init_func=init,frames=100, interval=100)    
ani.save("./model/{}.mp4".format(file_name), writer="ffmpeg",dpi=400)
