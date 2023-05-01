#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <random>

using namespace std;

long double pi = 3.14159265358979323846;

long double det(vector<vector<long double>> mat) {
    if (mat.size() == 0) {
        return 0;
    }
    
    if (mat.size() == 1) {
        return mat[0][0];
    }
    
    if (mat.size() == 2) {
        return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
    }
    
    long double res = 0;
    for (size_t i = 0; i < mat.size(); i++) {
        vector<vector<long double>> sub;
        for (size_t r = 1; r < mat.size(); r++) {
            vector<long double> sub_row;
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

vector<vector<long double>> matmul (vector<vector<long double>> A, vector<vector<long double>> B) {
    size_t A_rownum = A.size();
    size_t A_colnum = A[0].size();
    
    size_t B_rownum = B.size();
    size_t B_colnum = B[0].size();
    
    vector<vector<long double>> res;
    
    if (A_colnum != B_rownum) {
        printf("Matrix dimension does not match\n");
        return res;
    }
    
    for (size_t i = 0; i < A_rownum; i++) {
        vector<long double> res_row;
        for (size_t j = 0; j < B_colnum; j++) {
            long double x = 0;
            for (size_t k = 0; k < A_colnum; k++) {
                x += A[i][k] * B[k][j];
            }
            res_row.push_back(x);
        }
        res.push_back(res_row);
    }
    
    return res;
}

void vec_add_inplace (vector<long double> &x, vector<long double> y) {
    vector<long double> res;
    
    if (x.size() != y.size()) {
        printf("Vector length does not match to minus\n");
    }
    
    for (size_t i = 0; i < x.size(); i++) {
        x[i] += y[i];
    }
}

vector<long double> vec_minus (vector<long double> x, vector<long double> y) {
    vector<long double> res;
    
    if (x.size() != y.size()) {
        printf("Vector length does not match to minus\n");
        return res;
    }
    
    for (size_t i = 0; i < x.size(); i++) {
        res.push_back(x[i] - y[i]);
    }
    
    return res;
}

vector<long double> vec_mul (long double x, vector<long double> y) {
    vector<long double> res;
    
    for (size_t i = 0; i < y.size(); i++) {
        res.push_back(x * y[i]);
    }
    
    return res;
}

vector<long double> vec_div (long double x, vector<long double> y) {
    vector<long double> res;
    
    for (size_t i = 0; i < y.size(); i++) {
        res.push_back(y[i] / x);
    }
    
    return res;
}

vector<vector<long double>> transpose (vector<vector<long double>> X) {
    size_t rownum = X.size();
    size_t colnum = X[0].size();
    
    vector<vector<long double>> res;
    for (size_t i = 0; i < colnum; i++) {
        vector<long double> r;
        for (size_t j = 0; j < rownum; j++) {
            r.push_back(X[j][i]);
        }
        res.push_back(r);
    }
    
    return res;
}

vector<vector<long double>> matrix_inverse(vector<vector<long double>> a) {
    int i, j, k;
    int n = a.size();
    long double d;
    vector<vector<long double>> b;
    vector<vector<long double>> res;
    
    for (i = 0; i < n; i++) {
        vector<long double> temp;
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
        vector<long double> res_row;
        for (j = 0; j < n; j++) {
            res_row.push_back(b[i][j + n]);
        }
        res.push_back(res_row);
    }
    
    return res;
}

vector<vector<long double>> identity(int dim) {
    vector<vector<long double>> res;
    for (int i = 0; i < dim; i++) {
        vector<long double> row;
        for (int j = 0; j < dim; j++) {
            row.push_back(i == j ? 1.0 : 0.0);
        }
        res.push_back(row);
    }
    return res;
}

long double gaussian (vector<long double> x, vector<long double> mean, vector<vector<long double>> cov) {
    long double M = 2.0;
    long double scale = pow((2*pi), (-(long double)M/2.0)) * pow(det(cov), (-(long double)1/2.0));
    
    //printf("Scale partials %Lf %Lf\n", det(cov), pow(det(cov), (-(long double)1/2.0)));
    vector<long double> x_minus_mean = vec_minus(x, mean);
    vector<vector<long double>> x_mean {x_minus_mean};

    vector<vector<long double>> x_mean_T = transpose(x_mean);
    vector<vector<long double>> cov_inv = matrix_inverse(cov);
    vector<vector<long double>> oxxx = matmul(x_mean, cov_inv);

    vector<vector<long double>> mul_res = matmul(matmul(x_mean, cov_inv), x_mean_T);
    long double index = - (mul_res[0][0]) / 2;
    
    //printf("Scale is %Lf, index is %Lf, exp is %Lf\n", scale, index, exp(index));
    return scale * exp(index);
}

vector<vector<long double>> zeros (int x, int y) {
    vector<vector<long double>> res;
    for (int i = 0; i < x; i++) {
        vector<long double> temp(y, 0);
        res.push_back(temp);
    }
    
    return res;
}

vector<vector<long double>> inplace_mult (vector<vector<long double>> a, vector<vector<long double>> b) {
    int M = a.size();
    int N = a[0].size();
    vector<vector<long double>> res;
    for (int i = 0; i < M; i++) {
        vector<long double> row;
        for (int j = 0; j < N; j++) {
            row.push_back(a[i][j] * b[i][j]);
        }
        res.push_back(row);
    }
    return res;
}

vector<vector<long double>> scaler_mult (long double a, vector<vector<long double>> b) {
    int M = b.size();
    int N = b[0].size();
    vector<vector<long double>> res;
    for (int i = 0; i < M; i++) {
        vector<long double> row;
        for (int j = 0; j < N; j++) {
            row.push_back(a * b[i][j]);
        }
        res.push_back(row);
    }
    return res;
}

vector<vector<long double>> scaler_dev (long double a, vector<vector<long double>> b) {
    int M = b.size();
    int N = b[0].size();
    vector<vector<long double>> res;
    for (int i = 0; i < M; i++) {
        vector<long double> row;
        for (int j = 0; j < N; j++) {
            row.push_back(b[i][j] / a);
        }
        res.push_back(row);
    }
    return res;
}

void scaler_dev_inplace (long double a, vector<vector<long double>> &b) {
    int M = b.size();
    int N = b[0].size();

    for (int i = 0; i < M; i++) {
        vector<long double> row;
        for (int j = 0; j < N; j++) {
            b[i][j] /= a;
        }
    }

}

vector<vector<long double>> inplace_add (vector<vector<long double>> a, vector<vector<long double>> b) {

    if (a.size() != b.size() || a[0].size() != b[0].size()) {
        printf("Wrong dimension to add with A %ld * %ld, B %ld * %ld\n", a.size(), a[0].size(), b.size(), b[0].size());
    }
    
    int M = a.size();
    int N = a[0].size();
    vector<vector<long double>> res;
    for (int i = 0; i < M; i++) {
        vector<long double> row;
        for (int j = 0; j < N; j++) {
            row.push_back(a[i][j] + b[i][j]);
        }
        res.push_back(row);
    }
    return res;
}

vector<vector<long double>> inplace_subtract (vector<vector<long double>> a, vector<vector<long double>> b) {
    int M = a.size();
    int N = a[0].size();
    vector<vector<long double>> res;
    for (int i = 0; i < M; i++) {
        vector<long double> row;
        for (int j = 0; j < N; j++) {
            row.push_back(a[i][j] - b[i][j]);
        }
        res.push_back(row);
    }
    return res;
}