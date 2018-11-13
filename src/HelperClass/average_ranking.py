silhouette = ['AHC-average', 'PCA', 'kmeans', 'kmedoids', 'SC-kmeans', 'AHC-single', 'BIRCH', 'AP', 'DBSCAN', 'OPTICS', 'SC-eigen']
gamma = ['PCA', 'AP', 'AHC-average', 'kmedoids', 'kmeans', 'BIRCH', 'SC-kmeans', 'OPTICS', 'DBSCAN', 'SC-eigen', 'AHC-single']
dbindex = ['PCA', 'kmeans', 'AHC-single', 'kmedoids', 'AHC-average', 'BIRCH', 'SC-kmeans', 'SC-eigen', 'DBSCAN', 'OPTICS', 'AP']
validity = ['DBSCAN', 'AHC-single', 'PCA', 'SC-kmeans', 'BIRCH', 'AHC-average', 'kmedoids', 'AP', 'OPTICS', 'kmeans', 'SC-eigen']

sil_norm = ['d_R', 'd_S', 'd_P', 'd_E', 'd_M', 'd_F', 'd_H', 'd_G']
gamma_norm = ['d_H', 'd_E', 'd_F', 'd_G', 'd_M', 'd_S', 'd_R', 'd_P']
db_norm = ['d_G', 'd_M', 'd_S', 'd_H', 'd_P', 'd_E', 'd_F', 'd_R']
validity_norm = ['d_R', 'd_M', 'd_P', 'd_F', 'd_H', 'd_S', 'd_E', 'd_G']

average_ranking = {clustering:0 for clustering in silhouette}
average_norm = {norm:0 for norm in sil_norm}

for clustering in average_ranking.keys():
	for evaluation in [silhouette, gamma, dbindex, validity]:
		average_ranking[clustering]+=evaluation.index(clustering)
	average_ranking[clustering]=float(average_ranking[clustering])/float(4.0)

order = [100, None]
for e in average_ranking.keys():
	if order[0]>=average_ranking[e]:
		order[0]=average_ranking[e]
		order[1]=e

print('For streamlines...')
print('Best clustering is {}: {}'.format(order[1],order[0]))
print(average_ranking)

for norm in average_norm.keys():
	for evaluation in [sil_norm, gamma_norm, db_norm, validity_norm]:
		average_norm[norm]+=evaluation.index(norm)
	average_norm[norm]=float(average_norm[norm])/float(4.0)

order = [100, None]
for e in average_norm.keys():
	if order[0]>=average_norm[e]:
		order[0]=average_norm[e]
		order[1]=e

print('Best norm is {}: {}'.format(order[1],order[0]))
print(average_norm)


print('For pathlines...')
silhouette = ['AHC-average', 'AHC-single', 'PCA', 'SC-kmeans', 'BIRCH', 'kmeans', 'kmedoids', 'DBSCAN', 'AP', 'OPTICS', 'SC-eigen']
gamma = ['AHC-average', 'DBSCAN', 'PCA', 'AHC-single', 'BIRCH', 'OPTICS', 'SC-kmeans', 'kmeans', 'kmedoids', 'AP', 'SC-eigen']
dbindex = ['PCA', 'AHC-single', 'BIRCH', 'AHC-average', 'kmedoids', 'kmeans', 'SC-kmeans', 'OPTICS', 'DBSCAN', 'SC-eigen', 'AP']
validity = ['AHC-single', 'DBSCAN', 'PCA', 'AHC-average', 'AP', 'OPTICS', 'SC-eigen', 'SC-kmeans', 'kmedoids', 'BIRCH', 'kmeans']

sil_norm = ['d_R', 'd_H', 'd_G', 'd_E', 'd_M', 'd_F', 'd_S', 'd_P']
gamma_norm = ['d_E', 'd_M', 'd_H', 'd_G', 'd_F', 'd_R', 'd_S', 'd_P']
db_norm = ['d_G', 'd_E', 'd_M', 'd_H', 'd_F', 'd_S', 'd_R', 'd_P']
validity_norm = ['d_R', 'd_M', 'd_E', 'd_G', 'd_S', 'd_H', 'd_F', 'd_P']
average_ranking = {clustering:0 for clustering in silhouette}
average_norm = {norm:0 for norm in sil_norm}

for clustering in average_ranking.keys():
	for evaluation in [silhouette, gamma, dbindex, validity]:
		average_ranking[clustering]+=evaluation.index(clustering)
	average_ranking[clustering]=float(average_ranking[clustering])/float(4.0)

order = [100, None]
for e in average_ranking.keys():
	if order[0]>=average_ranking[e]:
		order[0]=average_ranking[e]
		order[1]=e

print('Best clustering is {}: {}'.format(order[1],order[0]))
print(average_ranking)

for norm in average_norm.keys():
	for evaluation in [sil_norm, gamma_norm, db_norm, validity_norm]:
		average_norm[norm]+=evaluation.index(norm)
	average_norm[norm]=float(average_norm[norm])/float(4.0)

order = [100, None]
for e in average_norm.keys():
	if order[0]>=average_norm[e]:
		order[0]=average_norm[e]
		order[1]=e

print('Best norm is {}: {}'.format(order[1],order[0]))
print(average_norm)
