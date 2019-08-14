# BIRCH Clustering Algorithm
BIRCH is a linear time scanning implemented by B+ tree to merge all the candidates within the distance threshold into one cluster. There are several comments on this BIRCH implementation
- The original implementation has hardcoded dimension size, and it **needs to be revised** according to the input data sets
- The final clusters are **spherical shape** because of squared distance threshold is utilized
- It does not guarantee the clustering result will be meaningful since the final clustering result heavily relies on the input distance threshold
- In order for fair comparison with those clustering techniques with an input of cluster number, we developed a **binary search** algorithm with input of cluster numbers to decide the possible distance threshold
	- For spatial similarity measures, it can almost guarantee the **distance threshold** can be found to generate **approximate required number of clusters**
	- For some similarity measures, it can never find such distance threshold, hence we set the *max search iteration* to be 10 
