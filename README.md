# FlowCurveClustering

This code folder provides unsupervised machine learning techniques with similarity measures for clustering integral curves in flow visualization community. It is the released source code for our TVCG paper, **Integral Curve Clustering and Simplification for Flow Visualization: A Comparative Evaluation**, by L. Shi, R. S. Laramee and G. Chen (to appear).

## Author Information

**Author**: Lieyu Shi, University of Houston

**Email**: shilieyu91@gmail.com


## Implemented Clustering Algorithms

- **k-means**, with option of [k-means++](https://en.wikipedia.org/wiki/K-means%2B%2B)

- **PCA** clustering, borrowed from [Streamline Variability Plots for Characterizing the Uncertainty in Vector Field Ensembles](https://ieeexplore.ieee.org/abstract/document/7192675) (TVCG 2016)
	- The original paper adopts [average-linkage AHC](https://en.wikipedia.org/wiki/Hierarchical_clustering) as clustering the lower-dimensional representation of streamlines, but in our experiments we find [k-means](https://en.wikipedia.org/wiki/K-means_clustering) works better
	- Additionally, due to high overload of AHC, k-means is better recommended to be used for PCA-based clustering

- [k-medoids](https://en.wikipedia.org/wiki/K-medoids)

- **DBSCAN**, a conventional [density-based clustering algorithm](https://en.wikipedia.org/wiki/DBSCAN)
	- Parameter tuning is pretty difficult, and in implementation we either rely on user to set the parameters, or directly set minPts to be 6, and eps as the average 6-th smallest distance of the neighboring streamlines
	- It also provides user input for minPts and eps in case user feel interested. However, it is not easy to find the so-called optimal pair of parameters for each flow data sets. 

- [OPTICS](https://en.wikipedia.org/wiki/OPTICS_algorithm)
	- There is more parameter tuning to get better clustering effect than DBSCAN
	- It is more difficult to be implemented, especially for determing the final number of clusters. The standard way is [reference paper](http://www.dbs.ifi.lmu.de/Publikationen/Papers/OPTICS.pdf), and we implemente a simple version of the cluster number determination

- [BIRCH](https://en.wikipedia.org/wiki/BIRCH) from the [BIRCH: A New Data Clustering Algorithm and Its Applications](https://link.springer.com/content/pdf/10.1023/A:1009783824328.pdf)
	- The implementation is directly borrowed from the [C++ github implementation](https://github.com/fedyura/birch-clustering-algorithm)
	- The clustering requires a distance value as input, which is pretty hard to handle especially for unknown flow data sets. Instead, I implemented a `binary search` to find the approximate distance threshold with the input of how many clusters are required to generate
		- For example, if the required number of clusters is 10, the program will decide which distance threshold is the right input for BIRCH such that the finalized clusters can be around 10
	- **Important Note** Since the dimension of tree is statically fixed, every time the user should try to replace the dimensions with that of the customized data sets, in 
		- src/Birch/main.cpp
		- src/Birch/ClusterAnalysis.h
		- src/Birch/CFTree.h. 
		
		For example, if the new data set has 1000 dimensions for each line, please **manually replace** 4824u with 1000u in ___src/Birch/main.cpp___, ___src/Birch/ClusterAnalysis.h___, ___src/Birch/CFTree.h___

- [Agglomerative hierarchical clustering](https://en.wikipedia.org/wiki/Hierarchical_clustering) (AHC)
	- Provide **single**, **average** and **complete** linkage option
		- **Why Ward's method is not provided?**

			Because Ward's method has not been used in flow visualization
	- Optimized by a [tree pre-processing technique](https://www.cs.cornell.edu/~kb/publications/IRT08.pdf)

- **Agglomerative hierarchical clustering** (AHC) with distance (ahc_dist)
	- The only difference is that the input is distance threshold, similar to BIRCH
	- It is definitely more **difficult** to determine the parameter than the input of number of clusters

- [Spectral Clustering](https://en.wikipedia.org/wiki/Spectral_clustering) with two post-processing techniques (SC)
	- k-means as proposed by [the tutorial](https://www.cs.cmu.edu/~aarti/Class/10701/readings/Luxburg06_TR.pdf) and [On Spectral Clustering: Analysis and an Algorithm](https://ai.stanford.edu/~ang/papers/nips01-spectral.pdf) by Ng et al. (NIPS 2001)
	- Eigenrotation minimization proposed by [self-tuning spectral clustering](https://papers.nips.cc/paper/2619-self-tuning-spectral-clustering.pdf)
		- The local scaling factor is set 5% as suggested by [Blood Flow Clustering and Applications in Virtual Stenting of Intracranial Aneurysms](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=6702500) (TVCG 2014)

- [Affinity propagation](https://en.wikipedia.org/wiki/Affinity_propagation) (AP)
	- Linux binary of AP can be obtained in [Frey Lab webpage](http://genes.toronto.edu/index.php?q=affinity%20propagation). We use OpenMP to implement the C++ version similar to the [github sample](https://github.com/nojima/affinity-propagation-sparse) and test on simple point-cloud data set and get the exactly same result
	- However, Tao. et al. [FlowString paper](https://ieeexplore.ieee.org/document/6787131) (TVCG 2014) used two-level hierarchical affinity propagation for streamline segment clustering. So we also implemented this hierarchical AP clustering
		- The initial value is set to the minimal similarity as the preference value


## Similarity Measures

There are just so many similarity measures, and even before the submission we already tried several similarity measures in the code (**as many as 17**!). However, in the TVCG submission, we only compared these following similarity measures. The **number** in the parenthesis marks the similarity measure option in the code.

- Euclidean distance d_E (**0**)

- Fraction norm d_F (**1**) with p=0.5 for L_P norm
	- It is claimed to have better performance in high-dimensional space than Euclidean distance (L_2)

- Geometric similarity measure d_G (**2**) from [my own short paper](http://www2.cs.uh.edu/~chengu/Publications/3DFlowVis/curveClustering.pdf)
	- It requires equal number of points between two integral lines
	- The implementation needs to handle the corner case of **zero** cosine values

- Accumulated rotation difference d_R (**4**)
	- It is a streamline attribute based similarity measure
	
- Mean-of-Closest-Point (MCP) d_M (**12**)
	- It is the **state-of-the-art** similarity measure for streamlines (possibly pathlines)

- Hausdorff distance d_H (**13**) from [Streamline Embedding for 3D Vector Field Exploration](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5753894) (TVCG 2012)
	- It can preserve topologic structure in streamlines

- Signature-based similarity d_S (**14**) from [Similarity Measures for Enhancing Interactive Streamline Seeding](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6231627)
	- Since our input of streamlines/pathlines only have `discrete curvatures` as signature, we only use `discrete curvatures` for [Chi-test](https://en.wikipedia.org/wiki/Chi-squared_test)
	- It is a linear combination of chi-test of streamline signature and MCP distance
	- Number of segments for computing chi-test of signature is an arguable issue. 
		- To avoid confusion, we select a fixed value for the experimented data sets since our focus is not discussing the similarity measure itself, but the comparisons of different clustering techniques with different similarity measures.
		- Optimal parameter tuning is tedious and endless, so I let it go

- Adapted Procrustes distance d_P (**15**) proposed by [A Vocabulary Approach to Partial Streamline Matching and Exploratory Flow Visualization](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7117453) (TVCG 2016)
	- Calculated [Procrustes distance](https://en.wikipedia.org/wiki/Procrustes_analysis) between every 7 points of streamlines
		- Implementation is difficult. I directly write C++ code based on procrustes.m in Matlab and verify the result
	- The original paper performs streamline segment clustering and the distance between two streamlines are decided by complicated strategy similar to string distance
		- To make it simple, I directly choose the **average of Proscrutes distance among all segments** for two streamlines
	- Let go of parameter tuning and optimal parameter pairs.

- Time-series MCP d_T (**17**) by [Exploration of Blood Flow Patterns in Cerebral Aneurysms during the Cardiac Cycle](https://www.sciencedirect.com/science/article/pii/S0097849318300128?via%3Dihub) (C && G 2018)
	- It is a similarity measure for pathlines considering time matching and overlapping, and claimed better than MCP
	- In our implementation, since our input pathlines have **exact time matching** for each time step, then the calculation is simplified a lot
		- Users can try the original d_T with time overlapping and mismatching

##### Note that we ignore all the parameter tuning issues and **only consider the most basic parameter pairs**. Parameter tuning is always a nightmare for designing similarity measures in flow visualization every body tries to avoid, so I guess why MCP is still regards the state-of-the-art similarity measure is simply due to that **it is parameter-free**!


## Before running?
1. Should adjust the BIN_SIZE in src/Common/Metric.cpp/Line 3 which is related to Chi-test distance computation
2. Could adjust k (size of compressed eigen-vector) in src/SpectralClustering/SpectralClustering.cpp/Line 522
3. Should adjust item_type<600u> in src/BIRCH/main.cpp, src/BIRCH/ClusterAnalysis.h, to maxDimension of input dataset

## How to run?

```
sudo chmod +x build.sh

./build.sh

cd Release

./kmeans dataSetName(must be put in ../dataset/) dimension(e.g.,3)
```