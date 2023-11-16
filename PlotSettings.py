import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd

#Global Plot Setting

sns.set_style('ticks') # darkgrid, white grid, dark, white and ticks
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=16)     # fontsize of the x and y labels
plt.rc('axes', linewidth=1.65 )   # width of the frame
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('font', size=12)          # controls default text sizes

#Color palette for the thesis

RegionAlpha = 0.5
MainColor1 = (29/255, 36/255, 164/255)
MainColor2 = (155/255, 8/255, 0/255)
MainColor3 = (212/255, 58/255, 0)
BackgroundColor1 = (0.467,0.137,0.184)
BackgroundColor2 = (253/255, 191/255, 92/255)
BackgroundColor3 = (29/255, 62/255, 64/255)
Gray1 = (0.337,0.352,0.360)
Gray2 = (0.694,0.698,0.690)

