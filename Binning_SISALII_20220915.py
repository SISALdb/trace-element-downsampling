# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 19:55:56 2022

@author: skiba
"""
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
import math

def consecutive(data, stepsize=0):
    return np.split(np.r_[:len(data)], np.where(np.diff(data) != stepsize)[0]+1) 

def consecutive_vls(data, stepsize=0):
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

#all depths need to be increasing! #no empty row in depth array when importing !

#import isotopes with regular depth scale -> leave blank where isotope measurements are missing -> use .csv format
#first column = depths in mm, second column = values)
#no header
col_names2=['depth', 'iso']
a = pd.read_csv('GejkarCave2017_isodepthwd18O_wmissing.csv', names=col_names2, header=None, sep = ',', encoding='utf-8')
a = a[a['depth'].notna()]#remove nans in the depth -> but there shouldn't be any in the imported file!; comment file if downsampling should fail when there are nans in the depth array
tmpl_datavalues = a['iso']
tmpl_data = a['depth']
tmpl_data = np.append(tmpl_data, tmpl_data[tmpl_data.index[-1]] + (tmpl_data[tmpl_data.index[-1]]-tmpl_data[tmpl_data.index[-2]]))#assumes that sample thickness for last sample (highest depth) is the same as for previous sample and adds the depth to which the drill bit then extended (because in the following sample thickness is determined by the previous samples depth and the next samples depth)

#import trace element data (i.e. higher resolved data); first column = depths in mm, second column = values)
#no header
b = np.loadtxt('GejkarCave2017_LAdepthwMgCa.txt', comments='#', encoding='utf-8')
datavalues = b[:,1]
datadepths = b[:,0]

#this programme assumes continuously drilled isotopes!
res = np.round(tmpl_data[1:(len(tmpl_data))] - tmpl_data[0:-1], 2) #calculates sample thickness for each sample 
#get to know the data
plt.hist(res)
plt.yscale('log')
plt.show()

#the data will be separated into different splits according to sampling resolution -> this ensures that in the following for loop sections with different sampling resolution /sampling thickness, hiati and irregularities will be treated accordingly
splits = consecutive(res)
splits_vls = consecutive_vls(res)

binnedvalues = []
abc = []
for h in range(0, len(splits)):
    if len(splits[h]) > 1:
        for i in range(0, len(splits[h])):
            if i < len(splits[h])-1:
                x = np.where(datadepths < tmpl_data[splits[h][i+1]])
                y = np.where(datadepths >= tmpl_data[splits[h][i]])
                z = np.intersect1d(x, y)
                abc.append(z)
            else: 
                end = tmpl_data[splits[h]][i] - tmpl_data[splits[h]][i-1] #assumes that sample thickness for last sample (highest depth) is the same as for previous sample and adds the depth to which the drill bit then extended (because in the following sample thickness is determined by the previous samples depth and the next samples depth)
                x = np.where(datadepths < tmpl_data[splits[h][i]]+end)
                y = np.where(datadepths >= tmpl_data[splits[h][i]])
                z = np.intersect1d(x, y)
                abc.append(z)
    else:
        if splits_vls[h][0] > np.mean(res)*4: #hiati defined as bigger than 4 times the mean sampling resolution -> this needs to be changed
                x = np.where(datadepths < (tmpl_data[splits[h][0]]+(tmpl_data[splits[h]][0] - tmpl_data[splits[h-1]][-1]))) #there will be a problem if the first iso sample thicknes is an outlier! 
                y = np.where(datadepths >= tmpl_data[splits[h][0]])
                z = np.intersect1d(x, y)
                abc.append(z)
        else: #this is just a sampling thickness outlier but has not been defined as hiatus
                x = np.where(datadepths < tmpl_data[splits[h+1]][0])
                y = np.where(datadepths >= tmpl_data[splits[h][0]])
                z = np.intersect1d(x, y)
                abc.append(z)
for i in range(0, len(abc)):    
    d = np.mean([datavalues[abc[i]]])
    binnedvalues.append(d)

#check again how the sampled were drilled -> check how many splits there are and if everything worked etc
plt.scatter(tmpl_data[0:-1], res)    
plt.show()

#plot and check result    
fig, ax = plt.subplots()
ax.plot(tmpl_data[0:-1], tmpl_datavalues,color="blue")
ax.set_ylabel('$\mathregular{δ^{18}O}$ [‰ VPDB]', fontsize=14, color="blue")
ax.set_xlabel('Depth [mm]', fontsize=14)
#ax.set_ylim(-7, -9)
ax.grid(False)
ax1 = ax.twinx()
ax1.spines["right"].set_position(("axes", 1))
ax1.yaxis.tick_right()
ax1.yaxis.set_ticks_position("right")
ax1.yaxis.set_label_position("right") 
#ax.plot(binneddata["age"], binneddata["value"], color='red',linestyle='dashed')
#o = 0
#ax1.plot(locals()["binneddata_part"+str(o)]["depth"], locals()["binneddata_part"+str(o)]["value"], color="red")
#o = 1
ax1.plot(tmpl_data[0:-1], binnedvalues, color="green")
ax1.plot(datadepths, datavalues, color="black",alpha=0.3)
ax1.set_ylabel('Mg/Ca [mmol/mol]', fontsize=14, color="red")
ax1.spines['right'].set_visible(True)
ax1.grid(False)
plt.show()
#ax1.set_xlim(200, -5)#, 200)
#ax1.set_ylim(0,1.8) #PCa
#ax1.set_ylim(0,0.015) #BaCa
#ax1.set_ylim(0,0.00015) #UCa

ax1.set_ylim(0,4)#,10) # MgCa
ax1.set_ylim(-0.01,0.1) #SrCa

#save
c = [tmpl_data[0:-1], binnedvalues] 
with open("HaozhuCave_HZZ27_2016-2018_binned_SrCa_2022_10_01.txt", "w") as file:
    for x in zip(*c):
        file.write("{0}\t{1}\n".format(*x))