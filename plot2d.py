from sys import stdin
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
from scipy.interpolate import griddata
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as mcolors
import sys 

#プローブ　試料　モード　曲率　高さ　x　y　長さ 間隔
args = sys.argv
probe_model = args[1]
surface_model = args[2]
scan_type = args[3]
probe_radius = args[4]
scan_height = args[5]
scan_x = args[6]
scan_y = args[7]
molecule_height = args[8]
molecule_distance =  args[9]
NUMBER = args[10]
args = args[1:]
FILE_NAME = '_'.join(args) + '.jpg'
DIR_NAME = "./imgs/" + surface_model + "_" + probe_model + "_" + probe_radius + "/" + scan_height +"/" + molecule_distance + "_" + molecule_height
os.makedirs(DIR_NAME, exist_ok=True)

N = 50
xy = []
for i in range(N):
    x = list(map(float, input().split()))
    xy.append(x)
xy = np.array(xy)
fig = plt.figure(figsize=(0.5,0.5))
plt.axis("off")
norm =  mcolors.TwoSlopeNorm(vmin=-1*abs(max(xy.min(),xy.max(),key=abs)),vcenter=0., vmax=abs(max(xy.min(),xy.max(),key=abs)))
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
plt.imshow(xy,norm=norm,cmap='gray')
plt.savefig(DIR_NAME + "/" + scan_x + "_" + scan_y + "_" + NUMBER + ".jpg")
