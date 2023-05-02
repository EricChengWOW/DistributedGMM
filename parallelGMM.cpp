#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <random>
#include "mpi.h"
#include "linearAlgebra.h"
#include "timing.h"

using namespace std;

int K = 4;                                       // number of clusters
int dim = 2;
vector<vector<long double>> mu_list;             // List of Mean Vectors
vector<vector<vector<long double>>> sigma_list;  // List of Covariance Matrixes
vector<long double> pi_list;                     // List of Pi Vector (Probability for each Gaussian Distribution)
vector<vector<long double>> dataset;             // Matrix of per node Dataset
vector<vector<long double>> whole_dataset;       // Matrix of the whole dataset

// Initialize the global parameters
void init() {

    for (int k = 0; k < K; k++) {
        vector<long double> tmp1;
        for (int i = 0; i < dim; i++) {
            tmp1.push_back(rand() % 20 - 10);
        }
        mu_list.push_back(tmp1);

        sigma_list.push_back(identity(dim));

        pi_list.push_back(1.0 / (long double)(K));
    }
}

// Each Iteration of GMM algorithm for updating mean, covariance and pi
void iterate(int pid, int nproc, int size) {
    int N = dataset.size();    // local dataset size
    int M = dataset[0].size(); // data dimension
    vector<vector<long double>> eta = zeros(K, N);
    vector<long double> local_eta_sum;
    MPI_Request send_request[nproc];
    
    // Compute Local eta value
    for (int i = 0; i < N; i++) {
        long double s = 0.0;
        for (int j = 0; j < K; j++) {
            s += gaussian(dataset[i], mu_list[j], sigma_list[j]) * pi_list[j];
        }
        
        //if (s == 0.0) s = 2.225E-307;
        
        for (int k = 0; k < K; k++) {
            eta[k][i] += gaussian(dataset[i], mu_list[k], sigma_list[k]) * pi_list[k] / s;
        }
    }
    
    // Compute Local eta sum
    for (int i = 0; i < K; i++) {
        long double sum = 0.0;
        for (int j = 0; j < N ; j++) {
            sum += eta[i][j];
        }
        
        local_eta_sum.push_back(sum);
    }
    
    // Master Node gather eta sum values and compute sums
    long double cluster_eta_sum[K];
    if (pid == 0) {
        vector<long double> node_eta [nproc];
        
        for (int i = 0; i < nproc; i++) {
            node_eta[i].resize(K);
        }
        
        for (int i = 1; i < nproc; i++) {
            MPI_Recv(node_eta[i].data(), K * sizeof(long double), MPI_BYTE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        
        // Compute eta sum
        for (int i = 0; i < K; i++) {
            cluster_eta_sum[i] = local_eta_sum[i];
            for (int j = 1; j < nproc; j++) {
                cluster_eta_sum[i] += node_eta[j][i];
            }
        }
    } else {
        MPI_Isend(local_eta_sum.data(), K * sizeof(long double), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &send_request[0]);
    }
    
    // Distribute eta sum
    MPI_Bcast(cluster_eta_sum, K * sizeof(long double), MPI_BYTE, 0, MPI_COMM_WORLD);
    
    // Calcute weights for mean calculation
    long double weights[K][N];
    for (int i = 0; i < K; i++) {
        for (int j = 1; j < N; j++) {
            weights[i][j] = eta[i][j] / cluster_eta_sum[i];
        }
    }
    
    // Calculate partial means
    vector<vector<long double>> local_means = zeros(K, M);
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < N; j++) {
            vec_add_inplace(local_means[i], vec_div(cluster_eta_sum[i], vec_mul(eta[i][j], dataset[j])));
        }
    }
    
    // Master Node gather mean data
    if (pid == 0) {
        vector<vector<long double>> node_mean [nproc];
        
        for (int i = 0; i < nproc; i++) {
            node_mean[i] = zeros(K, M);
        }
        
        for (int i = 1; i < nproc; i++) {
            for (int j = 0; j < K; j++) {
                MPI_Recv(node_mean[i][j].data(), M * sizeof(long double), MPI_BYTE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        
        // Compute real mean for next iteration
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < M; j++) {
                mu_list[i][j] = local_means[i][j];
            }
        }
        
        for (int i = 1; i < nproc; i++) {
            mu_list = inplace_add(mu_list, node_mean[i]);
        }
    } else {
        for (int i = 0; i < K; i++) {
            MPI_Isend(local_means[i].data(), M * sizeof(long double), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &send_request[0]);
        }
    }
    
    // Distribute mean
    for (int i = 0; i < K; i++) {
        MPI_Bcast(mu_list[i].data(), M * sizeof(long double), MPI_BYTE, 0, MPI_COMM_WORLD);
    }
    
    // Calculate Pi
    for (int k = 0; k < K; k++) {
        pi_list[k] = cluster_eta_sum[k] / size;
    }
    
    // Calculate local covariance
    for (int k = 0; k < K; k++) {
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
    }
    
    // Master gather covariance matrix
    if (pid == 0) {
        vector<vector<vector<long double>>> node_cov [nproc];
        
        for (int i = 0; i < nproc; i++) {
            for (int r = 0; r < K; r++) {
                node_cov[i].push_back(zeros(M, M));
            }
        }
        
        for (int i = 1; i < nproc; i++) {
            for (int j = 0; j < K; j++) {
                for (int r = 0; r < M; r++) {
                    MPI_Recv(node_cov[i][j][r].data(), M * sizeof(long double), MPI_BYTE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        
        for (int k = 0; k < K; k++) {
            for (int i = 1; i < nproc; i++) {
                sigma_list[k] = inplace_add(sigma_list[k], node_cov[i][k]);
            }
        }
        
        for (int i = 0; i < K; i++) {
            sigma_list[i] = scaler_dev(cluster_eta_sum[i], sigma_list[i]);
        }
    } else {
        for (int i = 0; i < K; i++) {
            for (int r = 0; r < M; r++) {
                MPI_Isend(sigma_list[i][r].data(), M * sizeof(long double), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &send_request[0]);
            }
        }
    }
    
    for (int i = 0; i < K; i++) {
        for (int r = 0; r < M; r++) {
            MPI_Bcast(sigma_list[i][r].data(), M * sizeof(long double), MPI_BYTE, 0, MPI_COMM_WORLD);
        }
    }
}

int main(int argc, char *argv[])
{
    //std::cout << "Parallel GMM Hello World!\n\n";
    
    int pid;
    int nproc;
  
    // Initialize MPI
    MPI_Init(&argc, &argv);
    // Get process rank
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    // Get total number of processes specificed at start of run
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  
    // Node Report
    if (pid == 0) printf("\n\nRunning %d nodes version\n", nproc);
    
    // Parsing Commands
    if (argc == 3) {
        if (pid == 0) {
            std::cout << "GMM dataset file: " << argv[1] << " \n";
            std::cout << "Target Cluster Count: " << argv[2] << " \n\n";
        }
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

        whole_dataset.push_back(row);
    }
    
    // Get local data
    size_t per_node = (whole_dataset.size() + nproc - 1) / nproc;
    for (size_t i = 0; i < per_node && i < whole_dataset.size(); i++) {
        dataset.push_back(whole_dataset[per_node * pid + i]);
    }
    
    printf("Node %d gets %ld entries\n", pid, dataset.size());
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (pid == 0) printf("--- Start Algorithm ---\n\n");
    //auto base_time = steady_clock::now();
    Timer totalSimulationTimer;
    
    // Initialize Parameters
    dim = dataset[0].size();
    init();
    
    // Communiate Node Sizes
    int N = dataset.size();
    int sizes[nproc];
    int total_size = N;
    MPI_Request send_request[nproc];
    
    for (int i = 0; i < nproc; i++) {
        if (i != pid) {
            MPI_Isend(&N, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &send_request[i]);
        }
    }
    
    for (int i = 0; i < nproc; i++) {
        if (i != pid) {
            MPI_Recv(&sizes[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_size += sizes[i];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    printf("Size Communicated %d : %d\n", pid, total_size);
    
    int epoch = 100;
    for (int i = 0; i < epoch; i++) {
        // Algorithm iteration
        iterate(pid, nproc, total_size);
        
        // Log intermediate results
        if (i % 10 == 0 && pid == 0) {
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
    
    if (pid == 0) {
        printf("total simulation time: %.6fs\n", totalSimulationTime);
    }
    
    if (pid == 0) {
        // Log final mean
        for (int j = 0; j < K; j++) {
            cout << "(" << mu_list[j][0] << "," << mu_list[j][1] << ")" << " ";
        }
        cout<<endl;
        
        // Store results to file
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
    }
    
    MPI_Finalize();
    
    return 0;
}
