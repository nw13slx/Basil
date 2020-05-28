import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
plt.switch_backend('agg')
from matplotlib.colors import LinearSegmentedColormap

def set_font():
    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'
    matplotlib.rcParams['font.size'] = 8

# https://matplotlib.org/3.1.0/gallery/color/named_colors.html

tabcolors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
             'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
