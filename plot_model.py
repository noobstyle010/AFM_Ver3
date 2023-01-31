import numpy as np
import glob 
from PIL import Image
import matplotlib.pyplot as plt
from scipy.stats.mstats import mquantiles
from mpl_toolkits.mplot3d import Axes3D
import sys

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
    12 : [C_SIGMA, C_EPSILON, "black"],
    16 : [O_SIGMA, O_EPSILON, "blue"],
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

for i in range(np):
    x, y, z, n = map(float, input().split())
    if z<=20 :
        n = int(n)
        ax.scatter(x, y, z, color=ATOMS[n][2],edgecolor='black',s=(fact*ATOMS[n][0])**2)
for i in range(ns):
    x, y, z, n = map(float, input().split())
    n = int(n)
    ax.scatter(x, y, z, color=ATOMS[n][2],edgecolor='black',s=(fact*ATOMS[n][0])**2)
plt.savefig("./model/{}.jpg".format(file_name))
