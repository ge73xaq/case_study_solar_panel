import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance
from matplotlib import cm
import imageio
import random

random.seed(123)

mean1 = [1, 2]
mean2 = [7, 6]
mean3 = [1, 7]
mean4 = [6, 1]
mean5 = [4, 4]

cov1 = [[1, 0], [0, 1]]
cov2 = [[1, 0], [0, 1]]
cov3 = [[1, 0.5], [0.5, 1]]
cov4 = [[1, 0.25], [0.25, 1]]
cov5 = [[1, 0], [0, 1]]

x1, y1 = np.random.multivariate_normal(mean1, cov1, 20).T
x2, y2 = np.random.multivariate_normal(mean2, cov2, 20).T
x3, y3 = np.random.multivariate_normal(mean3, cov3, 20).T
x4, y4 = np.random.multivariate_normal(mean4, cov4, 20).T
x5, y5 = np.random.multivariate_normal(mean5, cov5, 20).T

cluster1 = list(zip(x1, y1))
cluster2 = list(zip(x2, y2))
cluster3 = list(zip(x3, y3))
cluster4 = list(zip(x4, y4))
cluster5 = list(zip(x5, y5))

list_of_points = []
for i in cluster1 + cluster2 + cluster3 + cluster4 + cluster5:
    list_of_points.append(i)

centroid1 = (0, 0)
centroid2 = (6, 0)
centroid3 = (3, 3)
centroid4 = (0, 5)
centroid5 = (7, 2)

counter = 0
dict_old = 'test'
while counter < 100:
    dict = {(centroid1[0], centroid1[1]): [],
            (centroid2[0], centroid2[1]): [],
            (centroid3[0], centroid3[1]): [],
            (centroid4[0], centroid4[1]): [],
            (centroid5[0], centroid5[1]): []}
    for i in list_of_points:
        best_centroid = None
        min_dist = np.inf
        for j in list(dict.keys()):
            if distance.euclidean(i, j)  < min_dist:
                min_dist = distance.euclidean(i, j)
                dict[j].append(i)

    if dict_old == dict:
        break

    centroid1 = np.mean(dict[centroid1], axis=0)
    centroid2 = np.mean(dict[centroid2], axis=0)
    centroid3 = np.mean(dict[centroid3], axis=0)
    centroid4 = np.mean(dict[centroid4], axis=0)
    centroid5 = np.mean(dict[centroid5], axis=0)

    centroid1 = (centroid1[0], centroid1[1])
    centroid2 = (centroid2[0], centroid2[1])
    centroid3 = (centroid3[0], centroid3[1])
    centroid4 = (centroid4[0], centroid4[1])
    centroid5 = (centroid5[0], centroid5[1])

    viridis = cm.get_cmap('RdYlBu')
    col = 0
    for j in dict.keys():
        plt.scatter(j[0], j[1], color='black')
        for i in dict[j]:
            plt.scatter(i[0], i[1], color=viridis(col/5))
        col += 1

    plt.title(f'Iteration {counter}')
    plt.savefig(f'/home/lukas/Schreibtisch/plot{counter}.png')
    plt.show()
    dict_old = dict.copy()

    counter += 1

filenames = [f'/home/lukas/Schreibtisch/plot{i}.png' for i in range(counter)]
images = []
for i in filenames:
    images.append(imageio.imread(i))
kargs = {'duration': 2}
imageio.mimsave('/home/lukas/Schreibtisch/clustering.gif', images, **kargs)
