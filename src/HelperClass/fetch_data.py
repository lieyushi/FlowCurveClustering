from os import listdir
from os.path import isfile, join

# the script was used to read data from inside the folder obeying restrict order
# should have norm X to denote what norm it is using
# then would detect the value with respective tag
# summerize

def get_distance_limit(file_position):
	with open(file_position) as f:
		content = f.readlines()
	# you may also want to remove whitespace characters like `\n` at the end of each line
	content = [x.strip() for x in content] 
	distance_range = {}
	for x in content:
		if x!='':
			norm_pos = x.find('norm')
			pca_pos = x.find('PCA')

			if norm_pos==-1 and pca_pos==-1:
				continue
			if norm_pos!=-1:
				start_pos = norm_pos+5
				end_pos = start_pos
				while x[end_pos]!=' ' and x[end_pos]!=',':
					end_pos+=1
				norm_str = x[start_pos:end_pos]
			elif pca_pos!=-1:
				norm_str = 'PCA'

			range_pos = x.find('(max - min) is')
			start_ = range_pos+len('(max - min) is')
			while x[start_]==' ':
				start_+=1
			distance_range[norm_str] = float(x[start_:])
	return distance_range


# read optimal clustering and inside you've multiple clustering algorithm
def extract_evaluation_data(distance_range, data_folder):
	evaluation = {}
	norm_list = ['0','1','2','4','12','13','14','15']
	for d_folder in listdir(data_folder):
		readme = data_folder+'/'+d_folder+'/README'
		with open(readme) as r:
			content = r.readlines()
		norm_found = False
		norm=None
		evaluation[d_folder] = {}
		for val in norm_list:
			evaluation[d_folder][val] = {'silhouette':-10000.0, 'gamma':-10000.0, 'db index':-10000.0, 'validity':-10000.0, 'time':-10000.0}

		if d_folder=='kmeans':
			evaluation['PCA'] = {}
			for val in norm_list:
				evaluation['PCA'][val] = {'silhouette':-10000.0, 'gamma':-10000.0, 'db index':-10000.0, 'validity':-10000.0, 'time':-10000.0}

		norm_found = False
		for x in content:
			if x=='' or x=='\n':
				continue
			if norm_found is False:
				norm_pos = x.find('norm')
				Norm_pos = x.find('Norm:')
				pca_pos = x.find('PCA')
				if norm_pos==-1 and Norm_pos==-1 and pca_pos==-1:
					continue

				if norm_pos!=-1:
					end_pos = norm_pos+5
					while(x[end_pos]!=' ' and end_pos<=len(x)-1) and x[end_pos]!='\n':
						end_pos+=1
					norm = x[norm_pos+5:end_pos]
				elif Norm_pos!=-1:
					end_pos = Norm_pos+6
					while(x[end_pos]!=' ' and end_pos<=len(x)-1) and x[end_pos]!='\n':
						end_pos+=1
					norm = x[norm_pos+6:end_pos]
				elif pca_pos!=-1:
					norm = 'PCA'
				norm_found = True
			
			if norm_found is True and (norm in norm_list or norm=='PCA'):
				sil_pos = x.find('silhouette:')
				gamma_pos = x.find('statistic is:')
				dbindex_pos = x.find('DB index is:')
				validity_pos = x.find('measure is:')
				measurement_pos = x.find('measurement is:')

				if sil_pos!=-1:
					start_pos = sil_pos+len('silhouette:')+1
					while x[start_pos]==' ':
						start_pos+=1
					end_pos=start_pos
					while x[end_pos]!=',' and x[end_pos]!='\n' and end_pos<=len(x)-1 and x[end_pos]!=' ':
						end_pos+=1
					val_str = x[start_pos:end_pos]
					if val_str !='-nan' and val_str !='inf':
						if norm=='PCA':
							evaluation[norm]['0']['silhouette'] = float(x[start_pos:end_pos])
						else:
							evaluation[d_folder][norm]['silhouette'] = float(x[start_pos:end_pos])

				if gamma_pos!=-1:
					start_pos = gamma_pos+len('statistic is:')+1
					while x[start_pos]==' ':
						start_pos+=1
					end_pos=start_pos
					while x[end_pos]!=',' and x[end_pos]!='\n' and end_pos<=len(x)-1 and x[end_pos]!=' ':
						end_pos+=1
					val_str = x[start_pos:end_pos]
					if val_str !='-nan' and val_str !='inf':
						if norm=='PCA':
							evaluation[norm]['0']['gamma'] = float(x[start_pos:end_pos])
						else:
							evaluation[d_folder][norm]['gamma'] = float(x[start_pos:end_pos])

				if dbindex_pos!=-1:
					start_pos = dbindex_pos+len('DB index is:')+1
					while x[start_pos]==' ':
						start_pos+=1
					end_pos=start_pos
					while x[end_pos]!=',' and x[end_pos]!='\n' and end_pos<=len(x)-1 and x[end_pos]!=' ':
						end_pos+=1
					val_str = x[start_pos:end_pos]
					if val_str !='-nan' and val_str !='inf':
						if norm=='PCA':
							evaluation[norm]['0']['db index'] = float(x[start_pos:end_pos])
						else:
							evaluation[d_folder][norm]['db index'] = float(x[start_pos:end_pos])
					norm_found = False

				if validity_pos!=-1:
					start_pos = validity_pos+len('measure is:')+1
					while x[start_pos]==' ':
						start_pos+=1
					end_pos=start_pos
					while x[end_pos]!=',' and x[end_pos]!='\n' and end_pos<=len(x)-1 and x[end_pos]!=' ':
						end_pos+=1
					val_str = x[start_pos:end_pos]
					if val_str !='-nan' and val_str !='inf':
						if norm=='PCA':
							evaluation[norm]['0']['validity'] = float(x[start_pos:end_pos])/distance_range[norm]
						else:
							evaluation[d_folder][norm]['validity'] = float(x[start_pos:end_pos])/distance_range[norm]

				elif measurement_pos!=-1:
					start_pos = measurement_pos+len('measurement is:')+1
					while x[start_pos]==' ':
						start_pos+=1
					end_pos=start_pos
					while x[end_pos]!=',' and x[end_pos]!='\n' and end_pos<=len(x)-1 and x[end_pos]!=' ':
						end_pos+=1
					val_str = x[start_pos:end_pos]
					if val_str !='-nan' and val_str !='inf':
						if norm=='PCA':
							evaluation[norm]['0']['validity'] = float(x[start_pos:end_pos])/distance_range[norm]
						else:
							evaluation[d_folder][norm]['validity'] = float(x[start_pos:end_pos])/distance_range[norm]

				pca_time_tag = x.find('PCA+K_Means operation takes:')
				kmeans_time_tag = x.find('K-means on norm')
				kmedoid_time_tag = x.find('Direct K_Means operation time for norm')

				if pca_time_tag!=-1 or kmeans_time_tag!=-1 or kmedoid_time_tag!=-1:
					takes = x.find('takes:')
					if takes==-1:
						raise ValueError('Error for time search!')
					start_pos = takes+len('takes:')
					while x[start_pos]==' ':
						start_pos+=1
					end_pos=start_pos
					while x[end_pos]!=',' and x[end_pos]!='\n' and end_pos<=len(x)-1 and x[end_pos]!=' ' and x[end_pos]!='s':
						end_pos+=1
					val_str = x[start_pos:end_pos]
					if norm=='PCA':
						evaluation[norm]['0']['time'] = float(x[start_pos:end_pos])
					else:
						evaluation[d_folder][norm]['time'] = float(x[start_pos:end_pos])
				else:
					takes = x.find('takes:')
					if takes!=-1:
						start_pos = takes+len('takes:')
						while x[start_pos]==' ':
							start_pos+=1
						end_pos=start_pos
						while x[end_pos]!=',' and x[end_pos]!='\n' and end_pos<=len(x)-1 and x[end_pos]!=' ' and x[end_pos]!='s':
							end_pos+=1
						val_str = x[start_pos:end_pos]
						if evaluation[d_folder][norm]['time']<=-9999.0:
							evaluation[d_folder][norm]['time'] = float(x[start_pos:end_pos])
						else:
							evaluation[d_folder][norm]['time'] += float(x[start_pos:end_pos])
		
	return evaluation


# read AP and sc_eigen that only has one README file inside
def extract_single_readme(distance_range, data_folder):
	evaluation = {data_folder:{}}
	norm_list = ['0','1','2','4','12','13','14','15']
	for val in norm_list:
		evaluation[data_folder][val] = {'silhouette':-10000.0, 'gamma':-10000.0, 'db index':-10000.0, 'validity':-10000.0, 'time':-10000.0}
		
	readme = data_folder+'/README'
	with open(readme) as r:
		content = r.readlines()
	norm_found = False
	norm=None

	norm_found = False
	for x in content:
		if x=='' or x=='\n':
			continue
		if norm_found is False:
			norm_pos = x.find('norm')
			Norm_pos = x.find('Norm:')
			if norm_pos==-1 and Norm_pos==-1:
				continue

			if norm_pos!=-1:
				end_pos = norm_pos+5
				while(x[end_pos]!=' ' and end_pos<=len(x)-1) and x[end_pos]!='\n':
					end_pos+=1
				norm = x[norm_pos+5:end_pos]
			elif Norm_pos!=-1:
				end_pos = Norm_pos+6
				while(x[end_pos]!=' ' and end_pos<=len(x)-1) and x[end_pos]!='\n':
					end_pos+=1
				norm = x[norm_pos+6:end_pos]
			norm_found = True

		if norm_found is True and norm in norm_list:
			sil_pos = x.find('silhouette:')
			gamma_pos = x.find('statistic is:')
			dbindex_pos = x.find('DB index is:')
			validity_pos = x.find('measure is:')
			measurement_pos = x.find('measurement is:')

			if sil_pos!=-1:
				start_pos = sil_pos+len('silhouette:')+1
				while x[start_pos]==' ':
					start_pos+=1
				end_pos=start_pos
				while x[end_pos]!=',' and x[end_pos]!='\n' and end_pos<=len(x)-1 and x[end_pos]!=' ':
					end_pos+=1
				val_str = x[start_pos:end_pos]
				if val_str !='-nan' and val_str !='inf':
					evaluation[data_folder][norm]['silhouette'] = float(x[start_pos:end_pos])

			if gamma_pos!=-1:
				start_pos = gamma_pos+len('statistic is:')+1
				while x[start_pos]==' ':
					start_pos+=1
				end_pos=start_pos
				while x[end_pos]!=',' and x[end_pos]!='\n' and end_pos<=len(x)-1 and x[end_pos]!=' ':
					end_pos+=1
				val_str = x[start_pos:end_pos]
				if val_str !='-nan' and val_str !='inf':
					evaluation[data_folder][norm]['gamma'] = float(x[start_pos:end_pos])

			if dbindex_pos!=-1:
				start_pos = dbindex_pos+len('DB index is:')+1
				while x[start_pos]==' ':
					start_pos+=1
				end_pos=start_pos
				while x[end_pos]!=',' and x[end_pos]!='\n' and end_pos<=len(x)-1 and x[end_pos]!=' ':
					end_pos+=1
				val_str = x[start_pos:end_pos]
				if val_str !='-nan' and val_str !='inf':
					evaluation[data_folder][norm]['db index'] = float(x[start_pos:end_pos])
				norm_found = False

			if validity_pos!=-1:
				start_pos = validity_pos+len('measure is:')+1
				while x[start_pos]==' ':
					start_pos+=1
				end_pos=start_pos
				while x[end_pos]!=',' and x[end_pos]!='\n' and end_pos<=len(x)-1 and x[end_pos]!=' ':
					end_pos+=1
				val_str = x[start_pos:end_pos]
				if val_str !='-nan' and val_str !='inf':
					evaluation[data_folder][norm]['validity'] = float(x[start_pos:end_pos])/distance_range[norm]

			elif measurement_pos!=-1:
				start_pos = measurement_pos+len('measurement is:')+1
				while x[start_pos]==' ':
					start_pos+=1
				end_pos=start_pos
				while x[end_pos]!=',' and x[end_pos]!='\n' and end_pos<=len(x)-1 and x[end_pos]!=' ':
					end_pos+=1
				val_str = x[start_pos:end_pos]
				if val_str !='-nan' and val_str !='inf':
					evaluation[data_folder][norm]['validity'] = float(x[start_pos:end_pos])/distance_range[norm]

			takes = x.find('takes:')
			if takes!=-1:
				start_pos = takes+len('takes:')
				while x[start_pos]==' ':
					start_pos+=1
				end_pos=start_pos
				while x[end_pos]!=',' and x[end_pos]!='\n' and end_pos<=len(x)-1 and x[end_pos]!=' ' and x[end_pos]!='s':
					end_pos+=1
				val_str = x[start_pos:end_pos]
				if evaluation[data_folder][norm]['time']<=-9999.0:
					evaluation[data_folder][norm]['time'] = float(x[start_pos:end_pos])
				else:
					evaluation[data_folder][norm]['time'] += float(x[start_pos:end_pos])

	return evaluation


# read some data files like birch, dbscan, optics
def extract_norm_readme(distance_range, data_folder):
	evaluation = {data_folder:{}}
	norm_list = ['0','1','2','4','12','13','14','15']
	for val in norm_list:
		evaluation[data_folder][val] = {'silhouette':-10000.0, 'gamma':-10000.0, 'db index':-10000.0, 'validity':-10000.0, 'time':-10000.0}
	
	for norm in listdir(data_folder):
		readme = data_folder+'/'+norm+'/README'
		with open(readme) as r:
			content = r.readlines()

		for x in content:
			if x=='' or x=='\n':
				continue

			sil_pos = x.find('silhouette:')
			gamma_pos = x.find('statistic is:')
			dbindex_pos = x.find('DB index is:')
			validity_pos = x.find('measure is:')
			measurement_pos = x.find('measurement is:')

			if sil_pos!=-1:
				start_pos = sil_pos+len('silhouette:')+1
				while x[start_pos]==' ':
					start_pos+=1
				end_pos=start_pos
				while x[end_pos]!=',' and x[end_pos]!='\n' and end_pos<=len(x)-1 and x[end_pos]!=' ':
					end_pos+=1
				val_str = x[start_pos:end_pos]
				if val_str !='-nan' and val_str !='inf':
					evaluation[data_folder][norm]['silhouette'] = float(x[start_pos:end_pos])

			if gamma_pos!=-1:
				start_pos = gamma_pos+len('statistic is:')+1
				while x[start_pos]==' ':
					start_pos+=1
				end_pos=start_pos
				while x[end_pos]!=',' and x[end_pos]!='\n' and end_pos<=len(x)-1 and x[end_pos]!=' ':
					end_pos+=1
				val_str = x[start_pos:end_pos]
				if val_str !='-nan' and val_str !='inf':
					evaluation[data_folder][norm]['gamma'] = float(x[start_pos:end_pos])

			if dbindex_pos!=-1:
				start_pos = dbindex_pos+len('DB index is:')+1
				while x[start_pos]==' ':
					start_pos+=1
				end_pos=start_pos
				while x[end_pos]!=',' and x[end_pos]!='\n' and end_pos<=len(x)-1 and x[end_pos]!=' ':
					end_pos+=1
				val_str = x[start_pos:end_pos]
				if val_str !='-nan' and val_str !='inf':
					evaluation[data_folder][norm]['db index'] = float(x[start_pos:end_pos])
				norm_found = False

			if validity_pos!=-1:
				start_pos = validity_pos+len('measure is:')+1
				while x[start_pos]==' ':
					start_pos+=1
				end_pos=start_pos
				while x[end_pos]!=',' and x[end_pos]!='\n' and end_pos<=len(x)-1 and x[end_pos]!=' ':
					end_pos+=1
				val_str = x[start_pos:end_pos]
				if val_str !='-nan' and val_str !='inf':
					evaluation[data_folder][norm]['validity'] = float(x[start_pos:end_pos])/distance_range[norm]

			elif measurement_pos!=-1:
				start_pos = measurement_pos+len('measurement is:')+1
				while x[start_pos]==' ':
					start_pos+=1
				end_pos=start_pos
				while x[end_pos]!=',' and x[end_pos]!='\n' and end_pos<=len(x)-1 and x[end_pos]!=' ':
					end_pos+=1
				val_str = x[start_pos:end_pos]
				if val_str !='-nan' and val_str !='inf':
					evaluation[data_folder][norm]['validity'] = float(x[start_pos:end_pos])/distance_range[norm]

			takes = x.find('takes:')
			if takes!=-1:
				start_pos = takes+len('takes:')
				while x[start_pos]==' ':
					start_pos+=1
				end_pos=start_pos
				while x[end_pos]!=',' and x[end_pos]!='\n' and end_pos<=len(x)-1 and x[end_pos]!=' ' and x[end_pos]!='s':
					end_pos+=1
				val_str = x[start_pos:end_pos]
				if evaluation[data_folder][norm]['time']<=-9999.0:
					evaluation[data_folder][norm]['time'] = float(x[start_pos:end_pos])
				else:
					evaluation[data_folder][norm]['time'] += float(x[start_pos:end_pos])

	return evaluation


def get_average(lmethod_evaluation, sc_eigen_evaluation):
	average_evaluation = {}
	for clustering in lmethod_evaluation.keys():
		average_evaluation[clustering] = {}
		for norm in lmethod_evaluation[clustering].keys():
			average_evaluation[clustering][norm] = {}
			for eval_metric in lmethod_evaluation[clustering][norm].keys():
				average_val = 0.0
				effective = 0
				if lmethod_evaluation[clustering][norm][eval_metric]>=-9999.0:
					average_val+=lmethod_evaluation[clustering][norm][eval_metric]
					effective+=1
				if sc_eigen_evaluation[clustering][norm][eval_metric]>=-9999.0:
					average_val+=sc_eigen_evaluation[clustering][norm][eval_metric]
					effective+=1

				if effective==0:
					average_evaluation[clustering][norm][eval_metric] = -10000.0
				else:
					average_evaluation[clustering][norm][eval_metric]=average_val/effective
	return average_evaluation


def generate_text(evaluation_data, storage_name):
	storage = open(storage_name, 'w')
	clustering_algorithms = ['kmeans', 'kmedoids', 'AHC_single', 'AHC_average', 'birch', 'dbscan', 'optics', 'sc_kmeans', 'sc_eigen', 'AP', 'PCA']
	norm_order = ['0', '1', '2', '4', '12', '13', '14', '15']
	for clustering in clustering_algorithms:
		if clustering in evaluation_data.keys():
			first=[]
			second=[]
			for norm in norm_order:
				val = evaluation_data[clustering][norm]['silhouette']
				if val>=-9999.0:
					first.append(val)
				else:
					first.append('-')

				val = evaluation_data[clustering][norm]['gamma']
				if val>=-9999.0:
					first.append(val)
				else:
					first.append('-')

				val = evaluation_data[clustering][norm]['db index']
				if val>=-9999.0:
					second.append(val)
				else:
					second.append('-')

				val = evaluation_data[clustering][norm]['validity']
				if val>=-9999.0:
					second.append(val)
				else:
					second.append('-')
			for x in first:
				storage.write('%s ' % x)
			storage.write('\n')
			for x in second:
				storage.write('%s ' % x)
			storage.write('\n')


def generate_time(evaluation_data, storage_name):
	storage = open(storage_name, 'w')
	clustering_algorithms = ['kmeans', 'kmedoids', 'AHC_single', 'AHC_average', 'birch', 'dbscan', 'optics', 'sc_kmeans', 'sc_eigen', 'AP', 'PCA']
	norm_order = ['0', '1', '2', '4', '12', '13', '14', '15']
	for clustering in clustering_algorithms:
		if clustering in evaluation_data.keys():
			first=[]
			second=[]
			for norm in norm_order:
				val = evaluation_data[clustering][norm]['time']
				if val>=-9999.0:
					first.append(val)
				else:
					first.append('-')

			for x in first:
				storage.write('%s ' % x)
			storage.write('\n')
			for x in second:
				storage.write('%s ' % x)
			storage.write('\n')


def merge_two_dicts(first, second):
	result = first.copy()
	result.update(second)
	return result


def extract_full_data():
	distance_range=get_distance_limit('dist_range')

	print(distance_range)

	lmethod_evaluation = extract_evaluation_data(distance_range, 'optimal_clustering/lmethod')
	sc_eigen_evaluation = extract_evaluation_data(distance_range, 'optimal_clustering/sc_eigen_number')
	average_evaluation = get_average(lmethod_evaluation, sc_eigen_evaluation)
	print(average_evaluation['PCA']['0'])

	ap_evaluation = extract_single_readme(distance_range, 'AP')
	full_evaluation = merge_two_dicts(average_evaluation, ap_evaluation)

	sc_eigen_evaluation = extract_single_readme(distance_range, 'sc_eigen')
	full_evaluation = merge_two_dicts(full_evaluation, sc_eigen_evaluation)

	birch_evaluation = extract_norm_readme(distance_range, 'birch')
	full_evaluation = merge_two_dicts(full_evaluation, birch_evaluation)

	dbscan_evaluation = extract_norm_readme(distance_range, 'dbscan')
	full_evaluation = merge_two_dicts(full_evaluation, dbscan_evaluation)

	optics_evaluation = extract_norm_readme(distance_range, 'optics')
	full_evaluation = merge_two_dicts(full_evaluation, optics_evaluation)

	generate_text(full_evaluation, 'evaluation')

	generate_time(full_evaluation, 'time')


if __name__ == '__main__':
	extract_full_data()