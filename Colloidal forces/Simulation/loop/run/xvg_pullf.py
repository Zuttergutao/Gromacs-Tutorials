import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import sys


def distance(dat):
    with open(dat) as f:
        lines = f.readlines()
        labelName = [l for l in lines if l.startswith('@')]
        data = [[float(i) for i in l.split()] for l in lines if not (
            l.startswith('#') or l.startswith('@'))]

    title = re.findall(r'\"([\w\d\s\(\)\/]+)\"', labelName[0])[0]
    xlabel = re.findall(r'\"([\w\d\s\(\)\/]+)\"', labelName[1])[0]
    ylabel = re.findall(r'\"([\w\d\s\(\)\/]+)\"', labelName[2])[0]
    df = pd.DataFrame(data)
    plt.figure(figsize=(15, 8))
    plt.gca().spines["top"].set_linewidth(2)
    plt.gca().spines["bottom"].set_linewidth(2)
    plt.gca().spines["left"].set_linewidth(2)
    plt.gca().spines["right"].set_linewidth(2)
    plt.rcParams.update({"font.size": 24})
    df2 = df.drop([0], axis=1)
    plt.plot(df[0], df2)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlabel(xlabel, fontsize=24)
    plt.ylabel(ylabel, fontsize=24)
    plt.xlim(df[0].min(), df[0].max())
    plt.title(title, pad=12)
    plt.savefig('{}.png'.format(title),dpi=600)

if __name__ == '__main__':
    distance(sys.argv[1])
