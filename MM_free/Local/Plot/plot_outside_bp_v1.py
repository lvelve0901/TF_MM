
import numpy as np
import pandas as pd
import pylab as pl

barwidth = 0.2

free = np.array([66./409.,0./409.])
#bound = np.array([723./1736.,228./1736.])
tf = np.array([252./613.,69./613.])

#plot_list = [free,bound,tf]
#color_list = ['rosybrown','slateblue','limegreen']
#label_list = ['Free DNA','Protein DNA','TF DNA']

plot_list = [free,tf]
color_list = ['rosybrown','limegreen']
label_list = ['Free DNA','TF DNA']

fig = pl.figure(1,figsize=(10,5))
ax = fig.add_subplot(111)
ax.patch.set_facecolor('none')

for i in range(2):
    ax.bar(
        np.array(range(2))+barwidth*i-barwidth/2,
        plot_list[i],
        width = barwidth,
        color = color_list[i],
        linewidth=3,
        label=label_list[i]
    )

ax.set_xlim((0-len(plot_list)*barwidth/2-barwidth,1+len(plot_list)*barwidth/2+barwidth))
ax.set_xticks([0,1])
ax.set_xticklabels(["3xsigma","entire"],fontsize=20)
ax.set_ylim((0,0.8))
ax.set_ylabel("percentage",fontsize=20)

pl.legend(loc=1)
pl.savefig("./Plot/Survey_bp_params_outside_bar.pdf")

fig.tight_layout()
pl.show()

