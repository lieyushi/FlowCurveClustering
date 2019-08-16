## Spectral Clustering
Spectral clustering is a popular clustering technique that builds the normalized cut minimization for the input distance matrix. It includes two popular versions in flow visualization
- Spectral clustering (SC) with eigenrotation minimization (SC-eigen)
	- It can find the optimal number of clusters given the distance matrix and a preset bound k
	- It is very time consuming with complicated eigenrotation minimization inside the range
- k-means (SC k-means)
	- It finds the natural clusters with user input parameters after the generation of embedding space

Possible in the future we will implement a third version of spectral clustering, k-way normalized cut which has been found in flow visualization literature.
