# FlowCurveClustering
This provides an unsupervised machine learning techniques for clustering integral curves in flow visualization

------------------------------------------------
Implemented clustering techniques include

	K-means

	K-medoids

	DBSCAN (see link https://en.wikipedia.org/wiki/DBSCAN)

	OPTICS (see link https://en.wikipedia.org/wiki/OPTICS_algorithm)

	BIRCH (a fast hierarchical clustering)

	Agglomerative hierarchical clustering (ahc)(primitive version with one-pair-merge every time)

	Agglomerative hierarchical clustering with distance (ahc_dist) (updated version with pairs merged given distance of threshold)

	Spectral Clustering (1. k-means, 2. eigenvector rotation from self-tuning spectral clustering)

	The distance metric is provided from my paper (http://www2.cs.uh.edu/~chengu/Publications/3DFlowVis/curveClustering.pdf). Would incorporate more eixsting metrics.


-------------------------------------------------
How to run?

	sudo chmod +x build.sh

	./build.sh

	cd Release

	./kmeans dataSetName(must be put in ../dataset/) dimension(e.g.,3)

------------------------------------------------
Author information

	Lieyu Shi

	shilieyu91@gmail.com
	
	Computer Science @University of Houston, Houston, TX, USA
