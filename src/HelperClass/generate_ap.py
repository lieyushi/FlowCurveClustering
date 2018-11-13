streamline = 'birch/0/Crayfish_full.vtk'
ap_13 = 'AP/Crayfish_full.vtk'
result = 'ap_13.vtk'

with open(streamline) as f:
	content = f.readlines()
	# you may also want to remove whitespace characters like `\n` at the end of each line
	content = [x.strip() for x in content] 

with open(ap_13) as f:
	ap = f.readlines()
	ap = [x.strip() for x in ap] 


storage = open(result, 'w')
for x in content:
	if x.find('SCALARS group int 1')!=-1:
		break;
	if x!='':
		storage.write('%s\n' % x)

del content

start = False
count = 0;
for x in ap:
	if start is False:
		takes = x.find('SCALARS AP_norm13 int 1')
		start = (takes!=-1)
	if start is True and count<1528404:
		storage.write('%s\n' % x)
		count+=1
del ap




