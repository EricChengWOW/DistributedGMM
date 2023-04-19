#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>
#include <fstream>      // std::ifstream
#include <sstream>
#include <cmath>
#include <random>

using namespace std;

double pi = 3.14159265358979323846;
int K = 4; // number of clusters
vector<vector<double>> mu_list;
vector<vector<vector<double>>> sigma_list;
vector<double> pi_list;
vector<vector<double>> dataset;

template<class ForwardIt> ForwardIt max_element(ForwardIt first, ForwardIt last) {
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
                x += A[i][k] * B[k][j];
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

vector<vector<double>> identity(int dim) {
    vector<vector<double>> res;
    for (int i = 0; i < dim; i++) {
        vector<double> row;
        for (int j = 0; j < dim; j++) {
            row.push_back(i == j ? 1.0 : 0.0);
        }
        res.push_back(row);
    }
    return res;
}


double gaussian (vector<double> x, vector<double> mean, vector<vector<double>> cov) {

    // M = 2
    // scale = (2*np.pi)**(-M/2)*np.linalg.det(cov)**(-1/2)
    // return scale*np.exp(-(1/2)*(x-mu).T @ np.linalg.inv(cov) @ (x-mu))

    double M = 2;
    double scale = pow((2*pi), (-(double)M/2.0)) * pow(det(cov), (-(double)1/2.0));
    // cout<<pow(det(cov), (-1/2))<<endl<< endl;
    vector<double> x_minus_mean = vec_minus(x, mean);
    vector<vector<double>> x_mean {x_minus_mean};
    // cout << x_mean[0][0] << " " << x_mean[0][1] << endl << endl;
    vector<vector<double>> x_mean_T = transpose(x_mean);
    vector<vector<double>> cov_inv = matrix_inverse(cov);
    // cout << cov_inv[0][0] << " " << cov_inv[0][1] << endl << cov_inv[1][0] << " " << cov_inv[1][1] << endl << endl;
    vector<vector<double>> oxxx = matmul(x_mean, cov_inv);
    // cout<<oxxx[0][0]<<" "<<oxxx[0][1]<<endl;
    vector<vector<double>> mul_res = matmul(matmul(x_mean, cov_inv), x_mean_T);
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

vector<vector<double>> inplace_mult (vector<vector<double>> a, vector<vector<double>> b) {
    int M = a.size();
    int N = a[0].size();
    vector<vector<double>> res;
    for (int i = 0; i < M; i++) {
        vector<double> row;
        for (int j = 0; j < N; j++) {
            row.push_back(a[i][j] * b[i][j]);
        }
        res.push_back(row);
    }
    return res;
}

vector<vector<double>> scaler_mult (double a, vector<vector<double>> b) {
    int M = b.size();
    int N = b[0].size();
    vector<vector<double>> res;
    for (int i = 0; i < M; i++) {
        vector<double> row;
        for (int j = 0; j < N; j++) {
            row.push_back(a * b[i][j]);
        }
        res.push_back(row);
    }
    return res;
}

vector<vector<double>> scaler_dev (double a, vector<vector<double>> b) {
    int M = b.size();
    int N = b[0].size();
    vector<vector<double>> res;
    for (int i = 0; i < M; i++) {
        vector<double> row;
        for (int j = 0; j < N; j++) {
            row.push_back(b[i][j] / a);
        }
        res.push_back(row);
    }
    return res;
}

vector<vector<double>> inplace_add (vector<vector<double>> a, vector<vector<double>> b) {
    int M = a.size();
    int N = a[0].size();
    vector<vector<double>> res;
    for (int i = 0; i < M; i++) {
        vector<double> row;
        for (int j = 0; j < N; j++) {
            row.push_back(a[i][j] + b[i][j]);
        }
        res.push_back(row);
    }
    return res;
}

vector<vector<double>> inplace_subtract (vector<vector<double>> a, vector<vector<double>> b) {
    int M = a.size();
    int N = a[0].size();
    vector<vector<double>> res;
    for (int i = 0; i < M; i++) {
        vector<double> row;
        for (int j = 0; j < N; j++) {
            row.push_back(a[i][j] - b[i][j]);
        }
        res.push_back(row);
    }
    return res;
}

void init() {

    for (int k = 0; k < K; k++) {
        vector<double> tmp1;
        tmp1.push_back(1.0 * (k + 1));
        tmp1.push_back(1.0 * (k + 1));
        mu_list.push_back(tmp1);

        sigma_list.push_back(identity(2));

        pi_list.push_back(1.0 / double(K));
    }
}

void iterate() {
    int N = dataset.size();
    int M = dataset[0].size();
    vector<vector<double>> eta = zeros(K, N);
    for (int i = 0; i < N; i++) {
        double s = 0;
        for (int j = 0; j < K; j++) {
            s += gaussian(dataset[i], mu_list[j], sigma_list[j]) * pi_list[j];
            // if (i == 0 && j == 0) {
            //     cout << dataset[i][0] << " " << dataset[i][1] << endl;
            //     cout << mu_list[j][0] << " " << mu_list[j][1] << endl;
            //     cout << sigma_list[j][0][0] << " " << sigma_list[j][0][1] << endl;
            //     cout << sigma_list[j][1][0] << " " << sigma_list[j][1][1] << endl;
            //     cout << pi_list[j] << endl;
            // }
        }
        for (int k = 0; k < K; k++) {
            eta[k][i] += gaussian(dataset[i], mu_list[k], sigma_list[k]) * pi_list[k] / s;
        }
    }
    // for (int k = 0; k < K; k++) {
    //     for (int i = 0; i < N; i++) {
    //         cout << eta[k][i] << " ";
    //     }
    //     cout<<endl;
    // }
    for (int k = 0; k < K; k++) {
        double s = 0;
        for (int i = 0; i < N; i++) {
            s += eta[k][i];
        }
        pi_list[k] = s / N;
        vector<double> tmp1;
        for (int i = 0; i < M; i++) {
            tmp1.push_back(0.0);
        }
        for (int i = 0; i < N; i++) { 
            vector<double> tmp7;
            for (int j = 0; j < dataset[i].size(); j++) {
                tmp1[j] += eta[k][i] * dataset[i][j];
            }
        }
        for (int i = 0; i < tmp1.size(); i++) {
            tmp1[i] /= s;
        }
        mu_list[k] = tmp1;
        vector<vector<double>> tmp = zeros(M, M);
        for (int i = 0; i < N; i++) {
            vector<vector<double>> tmp2, tmp3, tmp4, tmp5, tmp6;
            tmp2.push_back(dataset[i]);
            tmp3.push_back(mu_list[k]);
            tmp4 = inplace_subtract(tmp2, tmp3);
            tmp5 = transpose(tmp4);
            tmp6 = matmul(tmp5, tmp4);
            sigma_list[k] = inplace_add(sigma_list[k], scaler_mult(eta[k][i], tmp6));
        }
        sigma_list[k] = scaler_dev(s, sigma_list[k]);
    }
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
    std::ifstream file(argv[1]);
    K = std::atoi(argv[2]);

    string line, field;

    // read data line by line
    while (std::getline(file, line))
    {
        vector<double> row;
        stringstream ss(line);

        // read field by field
        while (std::getline(ss, field, ','))
        {
            row.push_back(stod(field));
        }

        dataset.push_back(row);
    }
    
    printf("Recieve Dataset of Size %ld\n", dataset.size());

    // cout<<dataset[0].size() << endl;

    // Testing
    init();
    for (int i = 1; i < 101; i++) {
        iterate();
        if (i % 10 == 1) {
            for (int j = 0; j < K; j++) {
                cout << "(" << mu_list[j][0] << "," << mu_list[j][1] << ")" << " ";
            }
            cout<<endl;
        }
    }
    vector<double> x;
    // x.push_back(1.23123);
    // x.push_back(2.1233);
    // cout << x[0]<<" "<< x[1]<<endl;
    // vector<double> mean;
    // mean.push_back(1.2);
    // mean.push_back(0.223);
    // cout << mean[0]<<" "<< mean[1]<<endl;
    // vector<vector<double>> cov;
    // vector<double> tmp1, tmp2;
    // tmp1.push_back(1.56);
    // tmp1.push_back(0.902);
    // tmp2.push_back(0.123);
    // tmp2.push_back(0.4224);
    // cov.push_back(tmp1);
    // cov.push_back(tmp2);
    // cout << cov[0][0]<<" "<< cov[0][1]<<endl<< cov[1][0]<<" "<< cov[1][1]<<endl;
    // cout<<gaussian(x, mean, cov) << endl;
    return 0;
}
