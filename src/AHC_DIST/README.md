## AHC with distance threshold (**not frequently used in flow visualization**)
- It is very similar to [BIRCH](../Birch/README.md) which accepts a distance threshold and will merge all candidates into one group within this distance threshold. 
- It has not been applied to flow visualization and the intuition of implementing it is to compare the clustering result to that by Birch
- The implementation can be totally ignored
- The input distance threshold will be used to
	- If merged distance (calculated from the linkage type) is larger than the threshold, the hierarchical merging will terminate
	- If merged distance is below the threshold, the hierarchical merging will continue