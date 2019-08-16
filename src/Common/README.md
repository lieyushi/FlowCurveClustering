## Common Folder Description
It includes relevant functions 
- File I/O operations
- Sampling for streamlines/pathlines
- k-means initialization (from samples, from random coordinates, or k-means++)
- The hierarchical L-method for finding optimal number of clusters
- Different similarity measures for streamlines/pathlines
- The functions to calculate the clustering evaluation metrics, silhouette, the Gamma statics, DB index and normalized validity measurement

## Special notice

#### Distance Matrix
The distance matrix **distanceMatrix** is pre-stored as a 'float***' so that every time when calculating the similarity measure between two selected curves, the 'distanceMatrix' will be checked to be NULL or not. If 'distanceMatrix' is NULL, then the similarity measure function will be called otherwise the cached value is called.

#### MetricPreparation 
It is created before calculating the **MetricPreparation** due to the fact that for some similarity measures, e.g., the d_G (2), d_S(14) and d_P(15), either the segmentation on the streamlines/pathlines or the signature histograms should be calculated. In order to avoid repeated calculation of those signatures, we use a cache to pre-calculate the signatures for each line and store them for further pairwise distance value calculation.

