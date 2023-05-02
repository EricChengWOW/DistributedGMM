import matplotlib.pyplot as plt
import re
import numpy as np

def gaussian(x, mu, cov):
    # x and mu should be vectors in numpy, shape=(2,)
    # cov should be a matrix in numpy, shape=(2,2)
    M = 2
    scale = (2*np.pi)**(-M/2)*np.linalg.det(cov)**(-1/2)
    return scale*np.exp(-(1/2)*(x-mu).T @ np.linalg.inv(cov) @ (x-mu))

def plot_gaussian(mu, cov, x1_min=-100, x1_max=100, x2_min=-100, x2_max=100, color=None):
    # x and mu should be vectors in numpy, shape=(2,)
    # cov should be a matrix in numpy, shape=(2,2)

    x1_values = np.linspace(x1_min, x1_max, 101)
    x2_values = np.linspace(x2_min, x2_max, 101)

    x1_grid, x2_grid = np.meshgrid(x1_values,x2_values)

    M,N = x1_grid.shape
    y_grid = np.zeros((M,N))

    x = np.zeros((2,))

    for i in range(M):
        for j in range(N):
            x[0] = x1_grid[i,j]
            x[1] = x2_grid[i,j]

            y_grid[i,j] = gaussian(x, mu, cov)

    plt.contour(x1_grid, x2_grid, y_grid, colors=color)

if __name__ == "__main__":
    data = []
    xmin = 10000000
    xmax = -10000000
    ymin = 10000000
    ymax = -10000000
    
    with open("test_data.txt", "r") as log:
        lines = log.readlines()
        for line in lines:
            data_point = line.split(" ")
            data_point = data_point[0:len(data_point)-1]
            for i in range(len(data_point)):
                data_point[i] = float(data_point[i])
            # print(data_point)
            data.append(data_point)
            xmin = min(xmin, data_point[0])
            xmax = max(xmax, data_point[0])
            ymin = min(ymin, data_point[1])
            ymax = max(ymax, data_point[1])
    
    dim = 0
    mean = []
    cov = []
    to_mean = True
    with open("result.txt", "r") as log:
        lines = log.readlines()
        for line in lines:
            if line == "mean\n":
                continue
            if line == "\n":
                print("what")
                break
            if line == "cov\n":
                to_mean = False
                continue
            
            data_point = line.split(" ")
            data_point = data_point[0:len(data_point)-1]
            for i in range(len(data_point)):
                data_point[i] = float(data_point[i])
            
            if to_mean:
                dim = len(data_point)
                mean.append(data_point)
            else:
                cov_temp = []
                row = []
                for i in range(len(data_point)+1):
                    if i % dim == 0 and i != 0:
                        cov_temp.append(row)
                        row = []
                        if i == len(data_point):
                            break
                    row.append(data_point[i])
                cov.append(cov_temp)
                
    print(mean)
    print(cov)
    cluster_count = len(mean)
    per_cluster = len(data) // cluster_count

    for i in range(cluster_count):
        plt.scatter([d[0] for d in data[i * per_cluster : (i+1) * per_cluster]], [d[1] for d in data[i * per_cluster : (i+1) * per_cluster]], alpha=0.5)
    
    xRange = xmax - xmin
    yRange = ymax - ymin
    for i in range(len(mean)):
        plot_gaussian(mean[i], cov[i], xmin - xRange / 10, xmax + xRange / 10, ymin - yRange / 10, ymax + yRange / 10)
        
    plt.show()