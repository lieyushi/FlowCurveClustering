## Brief Descriptions for Folder Structures

As can be seen from the names of the folders, each folder represents the related calculations of one clustering technique. Some comments on the folders:

- **Common** [folder instruction](Common/README.md) includes: 
	- most of I/O operations
	- similarity measure calculation 
	- optimal cluster number search 
	- clustering evaluation metric calculation
	- a global pointer variable **called distanceMatrix**

- **HelperClass** [folder instruction](HelperClass/README.md) includes:
	- pathline resampling for blood flow that restrictly obeys the time matching criterion
		- direct repeating will be conducted during I/O processing of clustering algorithms
	- python scripts to extract all the analysis result data from the generated README by the experiments including
		- four evaluation metric values
		- calculation time
		- averaged rotation values for the cluster representative lines
		- final cluster information
	- C++ code to calculate the average and normalized evaluation metrics on all streamline and pathline data sets including
		- ranking-based sorting order of clustering techniques and similarity measures
		- normalized mapping of evaluation values into [0.1, 1.0]
		- calculate the average and standard deviation of evaluation values
		- output the latex tables for each and average evaluation results of streamline/pathline data sets
	- R code for ranking-based visualization
	- script to generate both analysis result and visualization plots
	- C++ code to read the limits of different similarity measures for each data set

- **Query** [folder instruction](Query/README.md) includes:
	- the streamline/pathline query results given a selected candidate 

- **ReadClustering** [folder instruction](ReadClustering/README.md) includes:
	- C++ code to read from the resulted .vtk file for clustering evaluation metric values just in case the evaluation values are out correctly output or saved in the clustering process

- **KMeans** [folder instruction](KMeans/README.md) includes:
	- PCA-based clustering
	- K-means with different similarity measures

- **Different clustering folders** includes
	- Clustering algorithm code for the given distance matrix
	- Clustering evaluation metric value calculation
	- Cluster representative calculation
	- In order to further accelerate the distance matrix caluation, for distance matrix of each similarity measure, if the **file** contains the distance matrix values in **../dataset/** (e.g., 0 file for Euclidean distance)
		- exists, it will **read the distance matrix from the local file**
		- does not exist, it will **calculate the distance matrix** from scratch and **store the distance values in the local file for further calculation**
