# DistributedGMM
Independent Project of Distributed Gaussian Mixture Model Clustering Algo

# Summary
	We are planning to implement the traditional Machine Learning clustering algorithm - Mixture of Gaussian in a distributed way. The potential parallelization in the problem includes data parallel for machine learning models, and SIMD processing for each distributed node.

# Background
  Clustering is a traditional important Machine Learning problem that focuses on grouping similar data points (samples, objects, etc) such that samples within the same cluster are more identical and with the similar features or characteristics.
  The Gaussian Mixture Model uses the EM concepts for a general 2-step algorithm. The expectation step represents a soft assignment for each data point to a cluster. The function is defined as a likelihood function for the current clustering status. 
  In this step, the contribution of each data point could be perfectly parallelized, but we would need a global synchronization to gather all intermediate parameters and create the whole likelihood function. 
  The maximization step involves the typical Machine Learning concept of MLE (Maximum Likelihood Estimation). Thus, given the function result of the E step, we will simply use the derivation of each GMM variable to gain the cluster parameters (mean, variance and mixture proportion for each sub-distribution) for the next iteration. 
From the distributed Perspective, this step either needs to be done sequentially with future parameter distribution to every node (parameter server analogy). Or we need all nodes to have the function, and compute the same next round simultaneously (ring communication analogy).

# Challenge
	Workload: 
Large scale dataset that can not be stored in a single computing resource
Complex computation required with huge amount of parameters in the EM algorithm

	Constraints:
In order to mimic the situation for multiple computing nodes or servers that each contains a subset of the data, we will use the message passing model which creates communication overhead
The challenge on the sequential portion inherent in the algorithm that could harm the speedup for the system.

# Resources
	We will start the project from scratch on the GHC machine with MPI to first create a valid model that behaves correctly. Then we will move to PSC with more available threads to measure the scalability and robustness of our algorithm.
	We do not need to refer to papers as the algorithm is already clear, but we do need to find benchmark datasets that we can use for clustering. In order to utilize the system on real world data, we might also need to encode or transform the data such that it could be processed by our system.

# Goals and Deliverables

	PLAN TO ACHIEVE:
Create the correct baseline GMM EM algorithm for clustering
Parallelize the algorithm with MPI to separate nodes/threads.
Test on different designs of communication models (parameter server, ring communication etc), and benchmark the comparison of performance

	HOPE TO ACHIEVE:
Add SIMD parallelization to each node for further utilization of the computing resource

	Demo:
The demo would contain the image illustration of how our algorithm performs and successfully cluster the dataset
The demo would presents the speedup graphs we have been creating for the labs in this course

	Analysis:
We hope to see the speedup for parallelizing the GMM clustering problem with high count of nodes
We hope to identify the bottleneck of the system for the potential non-linear speedup that we could have
We hope to measure the detailed work runtime, and use that to further confirm the Amdhal’s Law

# Platform
	We would like to use the GHC machine for development and PSC machine for final performance benchmarking since we want to test the scalability of our algorithm to a larger scale.
	We would develop the algorithm in C++ with access to MPI so that we can better mimic the behavior of multiple computing/storage nodes in the real world large scale situation.

# Schedule
4/3 - 4/9: Get more detailed information about all the computation process, the parallel algorithm, obtain the dataset, and starting to implement the baseline algorithm
4/10 - 4/16: Finish implementing the baseline GMM algorithm, come up with data visualization for how the algorithm clusters
4/17 - 4/23: Start to implement the distributed algorithm with parameter server
4/24 - 4/30: Finish implementing the parameter server version, start to implement the ring communication algorithm. *Try SIMD if time allows
5/1 - 5/5: Finish implementing both distributed versions, benchmark the performance, and create final analysis
