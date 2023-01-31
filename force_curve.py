import numpy as np
import glob 
from PIL import Image
import matplotlib.pyplot as plt
from scipy.stats.mstats import mquantiles
from mpl_toolkits.mplot3d import Axes3D

def turn_to_list(file_name):
    with open(file_name,encoding="utf-16") as f:
        forces = f.readlines()
    forces = [-1*float(force.rstrip('\n'))*6.94782*10 for force in forces]
    return forces

# Cur以下のものをグラフ化
save_path = "./force_curve/"
txt_path = "./cur/*"
# 2~7Å
Z = np.linspace(2,7,100)
Z = Z[20:-20]
files = glob.glob(txt_path)
for f in files:
    print(f[6:-4])
    F = turn_to_list(f)
    F = F[20:-20]
    plt.plot(Z,F)
    plt.hlines([np.min(F),np.min(F)],np.min(Z),np.max(Z),"black",linestyles="dashed")
    plt.xlabel("Z [Å]")
    plt.ylabel("Force [pN]")
    plt.xlim(np.min(Z),np.max(Z))
    plt.text(np.min(Z),np.min(F),"{:.1f}".format(np.min(F)),color='black')
    plt.savefig(save_path + f[6:-4]+".png")
    plt.cla()

