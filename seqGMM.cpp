#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>
#include <fstream>      // std::ifstream
#include <cmath>
#include <random>

using namespace std;
double pi = 3.14159265358979323846;

template<class ForwardIt>
ForwardIt max_element(ForwardIt first, ForwardIt last)
{
    if (first == last)
        return last;
 
    ForwardIt largest = first;
    ++first;
 
    for (; first != last; ++first)
        if (*largest < *first)
            largest = first;
 
    return largest;
}

double randn(double mean, double stddev) {
    static mt19937 gen{ random_device{}() };
    normal_distribution<> dist{ mean, stddev };
    return dist(gen);
}

double pdf(double x, double mean, double stddev) {
    return (1 / (stddev * sqrt(2 * pi))) * exp(-0.5 * pow((x - mean) / stddev, 2));
}

int argmax(const vector<double>& v) {
    return distance(v.begin(), max_element(v.begin(), v.end()));
}

vector<vector<double>> initialize_clusters(const vector<vector<double>>& X, int K) {
    vector<vector<double>> clusters(K, vector<double>(X[0].size()));
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < X[0].size(); j++) {
            double mean = X[rand() % X.size()][j];
            double stddev = 1;
            clusters[i][j] = randn(mean, stddev);
        }
    }
    return clusters;
}

vector<int> cluster_data(const vector<vector<double>>& X, int K) {
    vector<vector<double>> clusters = initialize_clusters(X, K);
    vector<int> labels(X.size());
    int iter = 1;
    
    while (true) {
        if (iter % 10000 == 0) printf("Iteration %d\n", iter);
        
        vector<vector<double>> responsibilities(X.size(), vector<double>(K));
        
        for (int i = 0; i < X.size(); i++) {
            for (int j = 0; j < K; j++) {
                double prior = 1.0 / K;
                double likelihood = pdf(X[i][0], clusters[j][0], clusters[j][1]);
                responsibilities[i][j] = prior * likelihood;
            }
            int label = argmax(responsibilities[i]);
            labels[i] = label;
        }
        
        vector<vector<double>> new_clusters(K, vector<double>(X[0].size()));
        vector<double> N(K);
        
        for (int i = 0; i < X.size(); i++) {
            int label = labels[i];
            for (int j = 0; j < X[0].size(); j++) {
                new_clusters[label][j] += X[i][j];
            }
            N[label] += 1;
        }
        
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < X[0].size(); j++) {
                new_clusters[i][j] /= N[i];
            }
        }
        
        double diff = 0;
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < X[0].size(); j++) {
                diff += abs(clusters[i][j] - new_clusters[i][j]);
                clusters[i][j] = new_clusters[i][j];
            }
        }
        
        if (diff < 1 || iter == 100000) {
            break;
        }
        
        iter++;
    }

    return labels;
}

int main(int argc, char *argv[])
{
    std::cout << "Sequential GMM Hello World!\n\n";
    
    // Parsing Commands
    if (argc == 3) {
        std::cout << "GMM dataset file: " << argv[1] << " \n";
        std::cout << "Target Cluster Count: " << argv[2] << " \n";
    } else {
        std::cout << "Not providing dataset or too many arguments\n";
        std::cout << "Format: ./seqGMM testfile.txt 5\n";
        return 1;
    }
    
    // Read dataset from file
    std::ifstream in(argv[1], std::ifstream::in);
    std::vector<std::vector<double>> dataset;
    int cluster_cnt = std::atoi(argv[2]);
    
    while(!in.eof()) {
      std::string datapoint;
      getline(in, datapoint);
      
      int start = 0;
      int end = 0;
      size_t size = datapoint.size();
      
      std::vector<double> datapoint_vec;
      for (size_t i = 0; i < size; i++) {
          if (datapoint.at(i) == ' ') {
              datapoint_vec.push_back(std::stod(datapoint.substr(start, i - start)));
              start = i+1;
          }
      }
      
      dataset.push_back(datapoint_vec);
      
    }
    dataset.erase(dataset.end());
    
    printf("Recieve Dataset of Size %ld\n", dataset.size());

    
    // Call the algorithm
    vector<int> labels = cluster_data(dataset, cluster_cnt);
    for (auto l : labels) {
        printf("%d\n", l);
    }
 
    return 0;
}
