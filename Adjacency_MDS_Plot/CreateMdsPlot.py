import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
#from matplotlib.collections import LineCollection

from sklearn import manifold
from sklearn.decomposition import PCA

dist_data_file = \
  "C:\\DePaulCoursework\\Research\\BryanExample\\2018715_21_44_53_Bicluster_Distance_Matrix.csv"
#dist_data_file = \
#  "C:\\DePaulCoursework\\Research\\BryanExample\\2018715_21_44_53_Bicluster_Adjacency_Low75.csv"
dist_df = pd.read_csv(dist_data_file, delimiter=",", index_col=0, header=0)
dist_np = np.array(dist_df)

#dist_np = dist_np.astype("float")
#ndx = dist_np == 0
#dist_np[ndx] = 0.000001
#dist_np = 1/dist_np

gene_list = list(dist_df.index)

seed = 20
mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=seed,
                   dissimilarity="precomputed", n_jobs=1)
pos = mds.fit(dist_np).embedding_

nmds = manifold.MDS(n_components=2, metric=False, max_iter=3000, eps=1e-12,
                    dissimilarity="precomputed", random_state=seed, n_jobs=1,
                    n_init=1)
npos = nmds.fit_transform(dist_np, init=pos)

# Rotate the data
clf = PCA(n_components=2)
pos = clf.fit_transform(pos)
npos = clf.fit_transform(npos)

fig = plt.figure(1)
ax = plt.axes([0., 0., 1., 1.])

s = 100
pos_x_data = pos[:, 0]
pos_y_data = pos[:, 1]
plt.scatter(pos_x_data, pos_y_data, color='turquoise', s=s, lw=0, label='MDS')
plt.legend(scatterpoints=1, loc='best', shadow=False)

for i, txt in enumerate(gene_list):
    ax.annotate(txt, (pos_x_data[i], pos_y_data[i]))

plt.show()
 
fig = plt.figure(1)
ax = plt.axes([0., 0., 1., 1.])

npos_x_data = npos[:, 0]
npos_y_data = npos[:, 1]
plt.scatter(npos_x_data, npos_y_data, color='darkorange', s=s, lw=0, label='NMDS')
plt.legend(scatterpoints=1, loc='best', shadow=False)

for i, txt in enumerate(gene_list):
    ax.annotate(txt, (npos_x_data[i], npos_y_data[i]))

plt.show()

#similarities = dist_np.max() / dist_np * 100
#similarities[np.isinf(similarities)] = 0

# Plot the edges
#start_idx, end_idx = np.where(pos)
# a sequence of (*line0*, *line1*, *line2*), where::
#            linen = (x0, y0), (x1, y1), ... (xm, ym)
#segments = [[X_true[i, :], X_true[j, :]]
#            for i in range(len(pos)) for j in range(len(pos))]
#values = np.abs(dist_np)
#lc = LineCollection(segments,
#                    zorder=0, cmap=plt.cm.Blues,
#                    norm=plt.Normalize(0, values.max()))
#lc.set_array(dist_np.flatten())
#lc.set_linewidths(np.full(len(segments), 0.5))
#ax.add_collection(lc)

