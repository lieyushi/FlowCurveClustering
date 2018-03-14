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

	Affinity propagation (very time comsuming, could only implement by O(n^3) but optimally could be O(kn^2 log{n}))

	The distance metric is provided from my paper (http://www2.cs.uh.edu/~chengu/Publications/3DFlowVis/curveClustering.pdf). Would incorporate more eixsting metrics. Finalized clustering metric number is 16.

-------------------------------------------------
Before running?

	1. Should adjust the BIN_SIZE in src/Common/Metric.cpp/Line 3 which is related to Chi-test distance computation
	2. Could adjust k (size of compressed eigen-vector) in src/SpectralClustering/SpectralClustering.cpp/Line 522
	3. Should adjust item_type<600u> in src/BIRCH/main.cpp, src/BIRCH/ClusterAnalysis.h, to maxDimension of input dataset

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
