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
N = 50
xy = []
for i in range(N):
    x = list(map(float, input().split()))
    xy.append(x)
np.array(xy)
plt.imshow(xy)
plt.show()