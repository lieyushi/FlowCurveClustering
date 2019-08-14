## Agglomerative Hierarachical Clustering (AHC)
The program includes basically two aspects of AHC
- AHC of three linkages (will generate cluster result information)
	- Single linkage
	- Complete linkage
	- Average linkage
- The hierarchical L method to find optimal number of clusters (only generate optimal cluster number)
	- It is a global search of knee point along the clusters

## Number of clusters as input
The program supports two kinds of input for number of clusters
- Direct input after the query information
	> Input cluster number among [0, 1000] for norm X: 
- Read the cluster numbers from a txt file
	- The txt file is called 'cluster_number' in the /dataset/ folder
	- The 'cluster number' has the following format

		0:10       // for similarity measure 0, the input of cluster number is 10
		1:10       // for similarity measure 1, the input of cluster number is 10
		2:10
		4:10
		12:10
		13:10
		14:10
		15:10
		17:10

		for better batch processing especially in our experiment when the code is automatically run on the server



