## Decription of Affinity Propagation
The implementation is an O(n^3) with OpenMP and we have tested the result on the point cloud data set and compre it to the [Frey Lab webpage](http://genes.toronto.edu/index.php?q=affinity%20propagation) linux binary version.

Two critical parameters are to be set
- 'Preference value' s(i,i)
	- The preference value in [affinity propagation](https://en.wikipedia.org/wiki/Affinity_propagation) and [Frey Lab webpage](http://genes.toronto.edu/index.php?q=affinity%20propagation) is set to be the **median** of negative squared Euclidean distance between points
	- However, in flow visualization, it is set to be the **minimal similarity value** among streamlines
- 'Relaxation factor' lambda
	- It controls the update rate and the default value is 0.5
- 'Max iteration'
	- Due to that distance matrix for the streamline data sets is often large size (>3000*3000), the default value is 20

## A two-level affinity propagation
Besides the conventional affinity propagation, the two-level affinity propagation is also included for user selection. 