import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
import os

def Histogram():
    total=pd.DataFrame()
    plt.figure(figsize=(35,8))
    for i in range(28):
        folders=os.getcwd()+"\prod_"+str(i)+"\md_pullx.xvg"
        with open(folders) as f:
            lines = f.readlines()
            labelName=[l for l in lines if l.startswith('@')]
            data=[[float(i) for i in l.split()] for l in lines if not (l.startswith('#') or l.startswith('@'))]
        df=pd.DataFrame(data)
        total.insert(loc=i,column=i,value=df[1])
        total[i].hist(stacked=True,bins=30,alpha=0.8,histtype="bar",edgecolor="black")
    plt.grid(None)
    plt.gca().spines["top"].set_linewidth(2)
    plt.gca().spines["bottom"].set_linewidth(2)
    plt.gca().spines["left"].set_linewidth(2)
    plt.gca().spines["right"].set_linewidth(2)
    plt.rcParams.update({"font.size":32})
    x_major_locator=MultipleLocator(0.1)
    plt.gca().xaxis.set_major_locator(x_major_locator)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.xlabel("Position (nm)",fontsize=32)
    plt.ylabel("Counts",fontsize=32)
    plt.title("Position Histogram",pad=24,fontsize=32)
    plt.savefig('{}.png'.format("Histogram"),dpi=600)

if __name__ == '__main__':
    Histogram()

