import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import lmfit
from itertools import cycle
import glob
from IPython.display import display
import seaborn as sns
sns.set()

size = 200
a2 = 1.0
a3 = 0.0
U = 1.0
df = pd.read_csv('data/dimer2_size' + str(size) + '_a2' + str(a2) + '_a3' + str(a3) + '_U' + str(U) + '.dat', sep="\t", header=None, index_col = 0)
#df = pd.read_csv('data/spin_size' + str(size) + '_U' + str(U) + '.dat', sep="\t", header=None, index_col = 0)
#df = pd.read_csv('data/entropy_mid_size' + str(size) + '.dat', sep="\t", header=None, index_col = 0)
x = df.index.to_numpy()
y = np.transpose(df[1].to_numpy())

x = x[0::1]
y = y[0::1]*(1)**(x)

pd.options.display.max_rows = None
pd.set_option('float_format', '{:.6f}'.format)

###

palette_dot = sns.color_palette('dark')
palette_line = sns.color_palette('dark')
palette_dot_iter = iter(palette_dot)
palette_line_iter = iter(palette_line)

#abc = next(palette_line_iter)
#abc = next(palette_dot_iter)
#abc = next(palette_line_iter)
#abc = next(palette_dot_iter)

with palette_line:
    for i in range(1):
        plt.plot(x, y, linewidth = 1, color=next(palette_line_iter), label='U = ' + str(U))

plt.xlim([x[0], x[-1]])
plt.xlabel("$r$")
plt.ylabel("$S$")
plt.legend(title = 'Entanglement Entropy(even)', loc='upper right')
#plt.gcf().subplots_adjust(left=0.2)
plt.show()
