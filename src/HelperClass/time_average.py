def get_average_time(time_list):
	data = []
	for each in time_list:
		data_set = []
		with open(each) as f:
			for line in f:
				if line!='' and line!='\n':
					line = line.strip()
					line = line.split()
					each_row = []
					for number in line:
						if number=='-':
							each_row.append(-10000.0)
						else:
							each_row.append(float(number))
					if each_row!=[]:
						data_set.append(each_row)
		if data_set!=[]:
			data.append(data_set)

	row = len(data[0])
	col = len(data[0][0])

	result = []
	for j in range(row):
		col_data = []
		for k in range(col):
			summation = 0 
			effective = 0
			for i in range(len(time_list)):
				if k<len(data[i][j]) and data[i][j][k]>=-9999.0:
					summation+=data[i][j][k]
					effective+=1
			if effective==0:
				summation = '-'
			else:
				summation=summation/effective

			col_data.append(summation)

		result.append(col_data)

	storage = open('average_time', 'w')
	for row in result:
		for x in row:
			storage.write('%s ' % x)
		storage.write('\n')
	

	
if __name__ == '__main__':
	streamlines = ['bernard_time', 'crayfish_time', 'cylinder_time', 'hurricane_time', 'solar_plume_time', 'tornado_time']
	pathlines = ['tub_pathlines_time', 'cylinder_pathlines_time', 'blood_flow_time']
	get_average_time(streamlines)