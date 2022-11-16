from sys import stdin
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
from scipy.interpolate import griddata

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as mcolors
import sys 

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
N = int(input())
for i in range(N):
    x, y, z, n = map(float, input().split())
    ax.scatter(x, y, z, color='purple')
plt.savefig("3D.jpg")