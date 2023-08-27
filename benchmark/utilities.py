
import matplotlib.transforms as mtransforms
from matplotlib.ticker import AutoMinorLocator
from matplotlib import pyplot as plt
import numpy as np

fontsize = 35
fontlegende = 35
font = {'family': 'sans', 'color':  'black', 'weight': 'normal', 'size': fontsize}
plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.serif": ["Palatino"]})

def complete_panel(ax, xlabel, ylabel, gray=[1, 1, 1], cancel_x=False, cancel_y=False, font=font, fontsize=fontsize, linewidth=2, tickwidth1=2.5, tickwidth2=2, legend=True, ncol=1, locator_x = 2, locator_y = 2, title=None):
    
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontdict=font)
        if cancel_x:
            ax.set_xticklabels([])
    else:
        ax.set_xticklabels([])

    if ylabel is not None:
        ax.set_ylabel(ylabel, fontdict=font)
        if cancel_y:
            ax.set_yticklabels([])  
    else:
        ax.set_yticklabels([])

    if title is not None:
        ax.set_title(title, fontdict=font)

    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    ax.yaxis.offsetText.set_fontsize(20)
    ax.minorticks_on()

    ax.tick_params('both', length=10, width=tickwidth1, which='major', direction='in', color=gray)
    ax.tick_params('both', length=6, width=tickwidth2, which='minor', direction='in', color=gray)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.label.set_color(gray)
    ax.yaxis.label.set_color(gray)

    ax.spines['left'].set_color(gray)
    ax.spines['top'].set_color(gray)
    ax.spines['bottom'].set_color(gray)
    ax.spines['right'].set_color(gray)

    ax.spines["top"].set_linewidth(linewidth)
    ax.spines["bottom"].set_linewidth(linewidth)
    ax.spines["left"].set_linewidth(linewidth)
    ax.spines["right"].set_linewidth(linewidth)

    ax.tick_params(axis='y', which='both', colors=gray)
    ax.tick_params(axis='x', which='both', colors=gray)

    minor_locator_x = AutoMinorLocator(locator_x)
    ax.xaxis.set_minor_locator(minor_locator_x)
    minor_locator_y = AutoMinorLocator(locator_y)
    ax.yaxis.set_minor_locator(minor_locator_y)

    if legend:
        leg = ax.legend(frameon=False, fontsize=fontlegende, 
                loc='best', handletextpad=0.5, ncol=ncol,
                handlelength = 0.86, borderpad = 0.3, 
                labelspacing=0.3)
        
        for text in leg.get_texts():
            text.set_color(gray)

myblue = np.array([20, 100, 255])/255
myred = np.array([255, 20, 20])/255
mygray = np.array([50, 50, 50])/255
myorange = np.array([255, 165, 0])/255
lightgray = [0.1, 0.1, 0.1]
darkgray = [0.9, 0.9, 0.9]