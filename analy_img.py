# 3d map 
import numpy as np
import glob 
from PIL import Image
import matplotlib.pyplot as plt
from scipy.stats.mstats import mquantiles
from mpl_toolkits.mplot3d import Axes3D

cb_min, cb_max = 0,0.01
cb_div = 10
interval_of_cf = np.linspace(cb_min, cb_max, cb_div+1)


Rs = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39, 40,41,42,43,44,45,46,47,48,49, 50,51,52,53,54,54,55,56,57,58,59, 60,61,62,63,64,65,66,67,68,69, 70]
probe_model = "FunctionalizedTip"
surface_model = "FourGold"
scan_mode = "Scan"
scan_height = "4"
scan_x = "0"
scan_y ="0"
molecule_length = 3
molecule_distances = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]

single_atom =  Image.open("./imgs/"+surface_model+"_"+"SingleGold"+"_-1/"+scan_height+"/-1_-1/"+scan_x+"_"+scan_y+"_"+"-1.jpg").convert(mode="L")
single_atom = np.array(single_atom)/255.0

max_dic = {}
min_dic = {}
ave_dic = {}

c=1
for a in molecule_distances:
    averages = []
    maxs = []
    mins = []
    for R in Rs:
        files = glob.glob("./imgs/"+surface_model+"_"+probe_model+"_"+str(R)+"/"+scan_height+"/"+str(a)+"_"+str(molecule_length)+"/"+scan_x+"_"+scan_y+"_"+"*")
        reference_img = np.zeros((50,50))
        indexes = []
        for file in files:
            img = Image.open(file).convert(mode="L") 
            img = np.array(img)/255.0
            index = np.sum(np.abs((img) - (single_atom))) / 2500
            indexes.append(index)
        averages.append(np.mean(indexes))
        maxs.append(np.max(indexes))
        mins.append(np.min(indexes))
    
        max_dic[(a,R)] = np.max(indexes)
        min_dic[(a,R)] = np.min(indexes)
        ave_dic[(a,R)] = np.mean(indexes)
##########################################################
X,Y = np.meshgrid(molecule_distances, Rs)
Z = [] 
for Ridx, R in enumerate(Rs):
    Z.append([])
    for aidx, a in enumerate(molecule_distances):
        Z[Ridx].append(max_dic[(a, R)])
Z = np.array(Z)
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
ax1.plot_surface(X, Y, Z,vmin=0,vmax=0.01, rstride=1, cstride=1, cmap='jet', linewidth=0.3) # 曲面のプロット。rstrideとcstrideはステップサイズ，cmapは彩色，linewidthは曲面のメッシュの線の太さ，をそれぞれ表す
ax1.set_zlim(0,0.05)
ax1.set_xlabel("a (Å)")
ax1.set_ylabel("R (Å)")
ax1.set_zlabel("MAE",rotation=90)
ax1.view_init(elev=60, azim=45)
plt.savefig("./FIGURE/3Dmax" + probe_model+surface_model+str(molecule_length)+"_"+scan_height+".jpg", bbox_inches="tight", pad_inches=0.1)
plt.cla()

plt.figure()
plt.contourf(X,Y,Z, interval_of_cf,cmap="jet",extend="both")
plt.xlabel("a (Å)")
plt.ylabel("R (Å)")
plt.savefig("./FIGURE/2Dmax" + probe_model+surface_model+str(molecule_length)+"_"+scan_height+".jpg")
plt.cla()
####################################################################
X,Y = np.meshgrid(molecule_distances, Rs)
Z = [] 
for Ridx, R in enumerate(Rs):
    Z.append([])
    for aidx, a in enumerate(molecule_distances):
        Z[Ridx].append(min_dic[(a, R)])
Z = np.array(Z)
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
ax1.plot_surface(X, Y, Z,vmin=0,vmax=0.01, rstride=1, cstride=1, cmap='jet', linewidth=0.3) # 曲面のプロット。rstrideとcstrideはステップサイズ，cmapは彩色，linewidthは曲面のメッシュの線の太さ，をそれぞれ表す
ax1.set_zlim(0,0.05)
ax1.set_xlabel("a (Å)")
ax1.set_ylabel("R (Å)")
ax1.view_init(elev=60, azim=45)
ax1.set_zlabel("MAE",rotation=90)
plt.savefig("./FIGURE/3Dmin" + probe_model+surface_model+str(molecule_length)+"_"+scan_height+".jpg", bbox_inches="tight", pad_inches=0.1)
plt.cla()

plt.figure()
plt.contourf(X,Y,Z, interval_of_cf,cmap="jet",extend="both")
plt.xlabel("a (Å)")
plt.ylabel("R (Å)")
plt.savefig("./FIGURE/2Dmin" + probe_model+surface_model+str(molecule_length)+"_"+scan_height+".jpg")
plt.cla()
####################################################################
X,Y = np.meshgrid(molecule_distances, Rs)
Z = [] 
for Ridx, R in enumerate(Rs):
    Z.append([])
    for aidx, a in enumerate(molecule_distances):
        Z[Ridx].append(ave_dic[(a, R)])
Z = np.array(Z)
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
ax1.plot_surface(X, Y, Z,vmin=0,vmax=0.01, rstride=1, cstride=1, cmap='jet', linewidth=0.3) # 曲面のプロット。rstrideとcstrideはステップサイズ，cmapは彩色，linewidthは曲面のメッシュの線の太さ，をそれぞれ表す
ax1.set_zlim(0,0.05)
ax1.set_xlabel("a (Å)")
ax1.set_ylabel("R (Å)")
ax1.view_init(elev=60, azim=45)
ax1.set_zlabel("MAE",rotation=90)
plt.savefig("./FIGURE/3Dave" + probe_model+surface_model+str(molecule_length)+"_"+scan_height+".jpg", bbox_inches="tight", pad_inches=0.1)
plt.cla()

plt.figure()
plt.contourf(X,Y,Z, interval_of_cf,cmap="jet",extend="both")
plt.xlabel("a (Å)")
plt.ylabel("R (Å)")
plt.savefig("./FIGURE/2Dave" + probe_model+surface_model+str(molecule_length)+"_"+scan_height+".jpg")
plt.cla()