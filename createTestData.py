import numpy as np

def create_cluster (dimension, num):
    mean = np.random.uniform(-5, 5, dimension)
    cov = np.random.normal(0, 2, (dimension, dimension))
    cov = np.eye(2)
    cov = np.dot(cov, cov.T)
    # cov = np.array([[6, -3], [-3, 3.5]])
    print(mean)
    
    return np.random.multivariate_normal(mean, cov, num)

if __name__ == "__main__":
    cluster_count = 5
    cluster_sizes = [100,100,100,100,100]
    dimension = 2
    
    data = []
    for i in range(cluster_count):
        data.append(create_cluster(dimension, cluster_sizes[i]))
    
    with open("test_data.txt", "w") as log:
        for c,cluster in enumerate(data):
            for i,d in enumerate(cluster):
                for feature in d:
                    log.write(str(feature))
                    log.write(" ")
                if (c != len(data)-1 or i != len(cluster)-1):
                    log.write("\n")
    
    total_count = sum(cluster_sizes)
    print("Created ", cluster_count, " clusters with total count of ", total_count)
            