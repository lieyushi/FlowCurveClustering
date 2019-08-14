## BIRCH Clustering Algorithm
BIRCH is a linear time scanning implemented by B+ tree to merge all the candidates within the distance threshold into one cluster. There are several comments on this BIRCH implementation
- The original implementation has hardcoded dimension size, and it **needs to be revised** according to the input data sets
- The final clusters are **spherical shape** because of squared distance threshold is utilized
- It does not guarantee the clustering result will be meaningful since the final clustering result heavily relies on the input distance threshold
- In order for fair comparison with those clustering techniques with an input of cluster number, we developed a **binary search** algorithm with input of cluster numbers to decide the possible distance threshold
	- For spatial similarity measures, it can almost guarantee the **distance threshold** can be found to generate **approximate required number of clusters**
	- For some similarity measures, it can never find such distance threshold, hence we set the *max search iteration* to be 10 


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