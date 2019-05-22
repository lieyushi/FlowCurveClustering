# FlowCurveClustering
This provides an unsupervised machine learning techniques for clustering integral curves in flow visualization community.
It is the released source code for one TVCG submission (to appear).

**Author**: Lieyu Shi, University of Houston,
**Email**: shilieyu91@gmail.com

1. **Implemented clustering algorithms**

	- K-means, see [k-means](https://en.wikipedia.org/wiki/K-means_clustering) 

	- PCA clustering, from [Streamline Variability Plots for Characterizing the Uncertainty in Vector Field Ensembles](https://ieeexplore.ieee.org/abstract/document/7192675)
		- The original paper adopts [average-linkage AHC](https://en.wikipedia.org/wiki/Hierarchical_clustering) as clustering for lower-dimensional space for streamlines, but in our experiments we find [k-means](https://en.wikipedia.org/wiki/K-means_clustering) works better
		- However, due to high overload of AHC, k-means is better to be used for PCA-based clustering

	- K-medoids, see [k-medoids](https://en.wikipedia.org/wiki/K-medoids)

	- DBSCAN, a conventional [density-based clustering algorithm](https://en.wikipedia.org/wiki/DBSCAN)
		- Parameter tuning is pretty difficult, and in implementation we either rely user to set the parameters, or directly set minPts to be 5, and eps as the average 5th smallest distance among all streamlines

	- OPTICS, see [OPTICS](https://en.wikipedia.org/wiki/OPTICS_algorithm)
		- More parameter tuning to get better clustering effect than DBSCAN
		- More difficult to be implemented, especially for determing the final number of clusters, see the [reference paper](http://www.dbs.ifi.lmu.de/Publikationen/Papers/OPTICS.pdf)
	
	- [BIRCH](https://en.wikipedia.org/wiki/BIRCH) from the [paper](https://link.springer.com/content/pdf/10.1023/A:1009783824328.pdf)
		- The implementation is directly borrowed from one [C++ github implementation](https://github.com/fedyura/birch-clustering-algorithm)
		- The clustering requires a distance value as input, which is pretty hard to handle especially for unknown data sets. Instead, I implemented a binary-search to find approximately cluster numbers if the input is how many clusters
		- **Important Note** Since the dimension of tree is statically fixed, every time the user should try to replace the dimensions with that of the customized data sets, in src/Birch/main.cpp, src/Birch/ClusterAnalysis.h, src/Birch/CFTree.h. 
		- For example, if the new data set has 1000 dimensions for each line, please manually replace 4824u with 1000u in src/Birch/main.cpp, src/Birch/ClusterAnalysis.h, src/Birch/CFTree.h

	- [Agglomerative hierarchical clustering](https://en.wikipedia.org/wiki/Hierarchical_clustering)
		- Provide **single**, **average** and **complete** linkage option
		- Optimized by a [tree pre-processing technique](https://www.cs.cornell.edu/~kb/publications/IRT08.pdf)

	- Agglomerative hierarchical clustering with distance (ahc_dist)
		- The only difference is that the input is distance threshold, similar to BIRCH

	- [Spectral Clustering](https://en.wikipedia.org/wiki/Spectral_clustering) with two post-processing techniques
		- K-means as proposed by [the tutorial](https://www.cs.cmu.edu/~aarti/Class/10701/readings/Luxburg06_TR.pdf) and [Andrew Y. Ng paper](https://ai.stanford.edu/~ang/papers/nips01-spectral.pdf)
		- Eigenrotation minimization proposed by [self-tuning spectral clustering](https://papers.nips.cc/paper/2619-self-tuning-spectral-clustering.pdf)

	- [Affinity propagation](https://en.wikipedia.org/wiki/Affinity_propagation)
		- Linux binary of AP can be obtained in [Frey Lab webpage](http://genes.toronto.edu/index.php?q=affinity%20propagation). We use OpenMP to implement the C++ version similar to the [github sample](https://github.com/nojima/affinity-propagation-sparse) and test on simple point-cloud data set and get the exactly same result
		- However, Tao. et al. [FlowString paper](https://ieeexplore.ieee.org/document/6787131) used two-level hierarchical affinity propagation for streamline segment clustering. So we finally implement this hierarchical AP clustering

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
