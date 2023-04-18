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

double det(vector<vector<double>> mat) {
    if (mat.size() == 0) {
        return 0;
    }
    
    if (mat.size() == 1) {
        return mat[0][0];
    }
    
    if (mat.size() == 2) {
        return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
    }
    
    double res = 0;
    for (size_t i = 0; i < mat.size(); i++) {
        vector<vector<double>> sub;
        for (size_t r = 1; r < mat.size(); r++) {
            vector<double> sub_row;
            for (size_t c = 0; c < mat.size(); c++) {
                if (c != i) sub_row.push_back(mat[r][c]);
            }
            sub.push_back(sub_row);
        }
        
        if (i % 2 ==0) {
            res += mat[0][i] * det(sub);
        } else {
            res -= mat[0][i] * det(sub);
        }
    }
    
    return res;
}

vector<vector<double>> matmul (vector<vector<double>> A, vector<vector<double>> B) {
    size_t A_rownum = A.size();
    size_t A_colnum = A[0].size();
    
    size_t B_rownum = B.size();
    size_t B_colnum = B[0].size();
    
    vector<vector<double>> res;
    
    if (A_colnum != B_rownum) {
        printf("Matrix dimension does not match\n");
        return res;
    }
    
    for (size_t i = 0; i < A_rownum; i++) {
        vector<double> res_row;
        for (size_t j = 0; j < B_colnum; j++) {
            double x = 0;
            for (size_t k = 0; k < A_colnum; k++) {
                x += A[i][k] + B[k][j];
            }
            res_row.push_back(x);
        }
        res.push_back(res_row);
    }
    
    return res;
}

vector<double> vec_minus (vector<double> x, vector<double> y) {
    vector<double> res;
    
    if (x.size() != y.size()) {
        printf("Vector length does not match to minus\n");
        return res;
    }
    
    for (size_t i = 0; i < x.size(); i++) {
        res.push_back(x[i] - y[i]);
    }
    
    return res;
}

vector<vector<double>> transpose (vector<vector<double>> X) {
    size_t rownum = X.size();
    size_t colnum = X[0].size();
    
    vector<vector<double>> res;
    for (size_t i = 0; i < colnum; i++) {
        vector<double> r;
        for (size_t j = 0; j < rownum; j++) {
            r.push_back(X[j][i]);
        }
        res.push_back(r);
    }
    
    return res;
}

vector<vector<double>> matrix_inverse(vector<vector<double>> a) {
    int i, j, k;
    int n = a.size();
    double d;
    vector<vector<double>> b;
    vector<vector<double>> res;
    
    for (i = 0; i < n; i++) {
        vector<double> temp;
        b.push_back(temp);
    }
    
    for (i = 0; i < n; i++) {
        for (j = 0; j < 2 * n; j++) {
            if (j < n) b[i].push_back(a[i][j]);
            else if (j == i + n) b[i].push_back(1);
            else b[i].push_back(0);
        }
    }
    
    for (i = 0; i < n; i++) {
        d = b[i][i];
        for (j = 0; j < 2 * n; j++) b[i][j] /= d;
        for (j = 0; j < n; j++) {
            if (i != j) {
                d = b[j][i];
                for (k = 0; k < 2 * n; k++) b[j][k] -= d * b[i][k];
            }
        }
    }
    
    for (i = 0; i < n; i++) {
        vector<double> res_row;
        for (j = 0; j < n; j++) {
            res_row.push_back(b[i][j + n]);
        }
        res.push_back(res_row);
    }
    
    return res;
}


double gaussian (vector<double> x, vector<double> mean, vector<vector<double>> cov) {
    double M = 2;
    double scale = pow((2*pi), (-M/2)) * pow(det(cov), (-1/2));
    
    vector<double> x_minus_mean = vec_minus(x, mean);
    vector<vector<double>> x_mean {x_minus_mean};
    vector<vector<double>> x_mean_T = transpose(x_mean);
    vector<vector<double>> cov_inv = matrix_inverse(cov);
    vector<vector<double>> mul_res = matmul(matmul(x_mean_T, cov_inv), x_mean);
     
    double index = - (mul_res[0][0]) / 2;
    return scale * exp(index);
}

vector<vector<double>> zeros (int x, int y) {
    vector<vector<double>> res;
    for (int i = 0; i < x; i++) {
        vector<double> temp(y, 0);
        res.push_back(temp);
    }
    
    return res;
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

    
    // Testing
    vector<vector<double>> a = {{1, 2, 1}, {3, 4, 1}, {1, 1, 1}};
    vector<vector<double>> inv = matrix_inverse(a);
    printf("det is %f\n", det(a));
    for (auto r : inv) {
        for (auto x : r) {
            printf("%f ", x);
        }
    }
    printf("\n");
 
    return 0;
}
