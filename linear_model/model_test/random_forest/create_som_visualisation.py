from mvpa2.suite import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle


# /usr/local/lib/python2.7/dist-packages/
sz = (10, 10)

p_agp = '/media/theo/Files/Lab/FHM/Complement/cross_valid_AGP_PRJEB/agp/'
p_prjeb = '/media/theo/Files/Lab/FHM/Complement/cross_valid_AGP_PRJEB/prjeb/'
p_inner = '/media/theo/Files/Lab/FHM/Complement/cross_valid_AGP_PRJEB/inner/'
p_small = '/media/theo/Files/Lab/FHM/Complement/cross_valid_AGP_PRJEB/inner_smallest/'

soms = []
for path in [p_small]:

    soms.append(SimpleSOMMapper(sz, 1000, learning_rate=0.05))

    with open(path + '/rf/conf_models.txt', 'r') as f:
        conf = f.read().replace('\r','').split('\n')
    conf.remove('')

    # conf = list(mods.index)[:-1]
    mods = pd.read_csv(path + '/log_errs/all_errs.txt', sep = '\t', header = 0, index_col= 0)
    conds = mods.iloc[-1,:]
    mods = mods.ix[conf,:].transpose()
    mods = mods.apply(pd.to_numeric, errors='ignore')

    mods_array = mods.as_matrix()
    soms[-1].train(mods_array)

som = soms[-1]
som._K = np.roll(soms[-1].K, 5, axis=0)
som._K = np.roll(som._K, -3, axis=1)
# plt.imshow(som2.K, origin='lower')
hc_heatmap  = np.zeros(sz)
cd_heatmap  = np.zeros(sz)

cd_mods = mods.loc[conds[conds == 'IBD'].index,:]
cd_hit = som.forward(cd_mods.as_matrix())
hc_mods = mods.loc[conds[conds == 'Healthy'].index,:]
hc_hit = som.forward(hc_mods.as_matrix())

for hit in hc_hit:
        hc_heatmap[hit[0],hit[1]] += 1
for hit in cd_hit:
        cd_heatmap[hit[0], hit[1]] += 1


# np.savetxt("foo.txt", bmu_heatmap, delimiter="\t")
plt.imshow(cd_heatmap/len(cd_hit), origin='lower')
plt.imshow(hc_heatmap/len(hc_hit), origin='lower')

slices = []
for var in range(som.K.shape[2]):
    slice = np.zeros(sz)
    for i in range(sz[0]):
        for j in range(sz[1]):
            slice[i,j] = som.K[i,j][var]
    # plt.imshow(slice, origin='lower')
    # plt.savefig('%i.png'%var)
    slices.append(slice)

###########
from mpl_toolkits.axes_grid1 import AxesGrid
from pylab import *

fig = plt.figure()

vals = slices

min_val = 100000
max_val = -100000

for s in slices:
    s_max = s.max()
    s_min = s.min()
    if min_val > s_min:
        min_val = s_min
    if max_val < s_max:
        max_val = s_max

grid = AxesGrid(fig, 111,
                nrows_ncols=(3,5),
                axes_pad=0.05,
                share_all=True,
                label_mode="L",
                cbar_location="right",
                cbar_mode="single",
                )

for val, ax in zip(vals,grid):
    im = ax.imshow(val, vmin=min_val, vmax=max_val)

grid.cbar_axes[0].colorbar(im)

for cax in grid.cbar_axes:
    cax.toggle_label(False)

# plt.axis('off')
# set_cmap('gray')
plt.show()

# savefig('inner_206vars.png')
som.K.tofile('12_12_15_smallest.out', sep='\t')
###########

# x = np.fromfile('12_12_15_smallest.out',sep='\t')
# x = x.reshape((12,12,15))

###########
# paint by dataset
agp_cols = list(pd.read_csv(p_agp + '/log_errs/all_errs.txt', sep ='\t').columns)[1:]
prjeb_cols = list(pd.read_csv(p_prjeb + '/log_errs/all_errs.txt', sep ='\t').columns)

agp_heatmap  = np.zeros(sz)
prjeb_heatmap  = np.zeros(sz)

ds = mods.iloc[-1,:]

agp_mods = mods.loc[[x for x in mods.index if x in agp_cols],:]
agp_hit = som.forward(agp_mods.as_matrix())
prjeb_mods = mods.loc[[x for x in mods.index if x in prjeb_cols],:]
prjeb_hit = som.forward(prjeb_mods.as_matrix())

for hit in agp_hit:
    agp_heatmap[hit[0],hit[1]] += 1
for hit in prjeb_hit:
    prjeb_heatmap[hit[0], hit[1]] += 1


# np.savetxt("foo.txt", bmu_heatmap, delimiter="\t")
plt.imshow(agp_heatmap/len(agp_hit), origin='lower')
plt.imshow(prjeb_heatmap/len(prjeb_hit), origin='lower')
