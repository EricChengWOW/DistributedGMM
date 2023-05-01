#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <random>
#include "linearAlgebra.h"
#include "timing.h"

using namespace std;

int K = 4; // number of clusters
vector<vector<long double>> mu_list;
vector<vector<vector<long double>>> sigma_list;
vector<long double> pi_list;
vector<vector<long double>> dataset;

void init() {

    for (int k = 0; k < K; k++) {
        vector<long double> tmp1;
        tmp1.push_back(1.0 * ((long double)k + 1) - (long double)K / 2.0);
        tmp1.push_back(1.0 * ((long double)k + 1) - (long double)K / 2.0);
        mu_list.push_back(tmp1);

        sigma_list.push_back(identity(2));

        pi_list.push_back(1.0 / (long double)(K));
    }
}

void iterate() {
    int N = dataset.size();
    int M = dataset[0].size();
    vector<vector<long double>> eta = zeros(K, N);
    
    for (int i = 0; i < N; i++) {
        long double s = 0.0;
        //long double s = 0.00001;
        for (int j = 0; j < K; j++) {
            s += gaussian(dataset[i], mu_list[j], sigma_list[j]) * pi_list[j];
        }
        
        //printf("%d %Lf\n", i, s);
        
        for (int k = 0; k < K; k++) {
            eta[k][i] += gaussian(dataset[i], mu_list[k], sigma_list[k]) * pi_list[k] / s;
            //printf("%Lf * %Lf / %Lf = %Lf\n", gaussian(dataset[i], mu_list[k], sigma_list[k]), pi_list[k], s, eta[k][i]);
        }
    }
    
    /*for (int k = 0; k < K; k++) {
        for (int i = 0; i < N; i++) {
             cout << eta[k][i] << " ";
        }
        cout<<endl;
    }*/
    
    for (int k = 0; k < K; k++) {
        // Compute the new Pi values
        long double s = 0.0;
        for (int i = 0; i < N; i++) {
            s += eta[k][i];
        }
        pi_list[k] = s / N;
        
        // Compute the new Mean values
        vector<long double> tmp1;
        for (int i = 0; i < M; i++) {
            tmp1.push_back(0.0);
        }
        
        for (int i = 0; i < N; i++) { 
            for (int j = 0; j < dataset[i].size(); j++) {
                tmp1[j] += eta[k][i] * dataset[i][j];
            }
        }
        
        for (int i = 0; i < tmp1.size(); i++) {
            tmp1[i] /= s;
        }
        mu_list[k] = tmp1;
        
        // Compute the new Sigma values
        sigma_list[k] = zeros(M, M);
        
        for (int i = 0; i < N; i++) {
            vector<vector<long double>> tmp2, tmp3, tmp4, tmp5, tmp6;
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
    std::cout << "\n\nSequential GMM Hello World!\n\n";
    
    // Parsing Commands
    if (argc == 3) {
        std::cout << "GMM dataset file: " << argv[1] << " \n";
        std::cout << "Target Cluster Count: " << argv[2] << " \n\n";
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
        vector<long double> row;
        stringstream ss(line);

        // read field by field
        while (std::getline(ss, field, ' '))
        {
            row.push_back(stod(field));
        }

        dataset.push_back(row);
    }
    
    printf("Recieve Dataset of Size %ld\n\n", dataset.size());

    for (int i = 0; i < 5; i++) {
        for (auto f : dataset[i]) {
            printf("%Lf ", f);
        }
        
        printf("\n\n");
    }
    
    Timer totalSimulationTimer;
    // Testing
    init();
    int epoch = 100;
    for (int i = 0; i < epoch; i++) {
        iterate();
        if (i % 10 == 0) {
            cout << "Mean at Iteration " << i << " ";
            for (int j = 0; j < K; j++) {
                cout << "(" << mu_list[j][0] << "," << mu_list[j][1] << ")" << " ";
            }
            cout<<endl;
            
            /*cout << "Pi at Iteration " << i << " ";
            for (int j = 0; j < K; j++) {
                cout << pi_list[j] << " ";
            }
            cout<<endl;
            
            cout << "Sigma at Iteration " << i << " ";
            for (int j = 0; j < K; j++) {
                cout << "(" << sigma_list[j][0][0] << "," << sigma_list[j][0][1] << ")" << " " << "(" << sigma_list[j][1][0] << "," << sigma_list[j][1][1] << ")";
            }
            cout<<endl<<endl<<endl;*/
        }
    }
    
    double totalSimulationTime = totalSimulationTimer.elapsed();
    
    printf("total simulation time: %.6fs\n", totalSimulationTime);
    
    for (int j = 0; j < K; j++) {
        cout << "(" << mu_list[j][0] << "," << mu_list[j][1] << ")" << " ";
    }
    cout<<endl;
    
    std::ofstream ofs ("result.txt", std::ofstream::out);

    ofs << "mean\n";
    
    for (auto mu : mu_list) {
        for (auto x : mu) {
            ofs << x << " ";
        }
        ofs << "\n";
    }
    
    ofs << "cov\n";
    
    for (auto sigma : sigma_list) {
        for (auto row : sigma) {
            for (auto x : row) {
                ofs << x << " ";
            }
        }
        ofs << "\n";
    }
    
    ofs.close();
    
    return 0;
}
