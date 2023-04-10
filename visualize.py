import matplotlib.pyplot as plt

if __name__ == "__main__":
    data = []
    with open("test_data.txt", "r") as log:
        lines = log.readlines()
        for line in lines:
            data_point = line.split(" ")
            data_point = data_point[0:len(data_point)-1]
            for i in range(len(data_point)):
                data_point[i] = float(data_point[i])
            # print(data_point)
            data.append(data_point)
    
    plt.scatter([d[0] for d in data], [d[1] for d in data], alpha=0.5)
    plt.show()