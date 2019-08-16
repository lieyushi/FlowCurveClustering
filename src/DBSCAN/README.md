## DBSCAN Description

#### Parameter selection
Two critical parameters for DBSCAN clustering
- The **minPts** which describes *how many neighbor candidates needed to create a core point*
- The **radius** which describes *how large the searched area can be to define a core point*
which are, however, pretty difficult for parameter tuning varying on different data sets. 

#### Our implementation details
We only use one parameter, **minPts**, and the other parameter, **radius**, is set to be the minPts-th smallest distance to the candidate line from all its neighbors.

**minPts** can be totally user defined, or to be default, set by 6 in all the data sets and it performs pretty well to generate around 10-100 clusters for fair comparisons of clustering combinations in our paper.

From this perspective, DBSCAN is not suitable for scientific data visualization depsite it has benefit of lower overhead and resource requirement for calculation.