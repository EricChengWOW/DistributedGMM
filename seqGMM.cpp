#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>
#include <fstream>      // std::ifstream

int main(int argc, char *argv[])
{
    std::cout << "Sequential GMM Hello World!\n\n";
    char *dataset_file;
    
    if (argc == 3) {
        std::cout << "GMM dataset file: " << argv[1] << " \n";
        std::cout << "Target Cluster Count: " << argv[2] << " \n";
    } else {
        std::cout << "Not providing dataset or too many arguments\n";
        std::cout << "Format: ./seqGMM testfile.txt 5\n";
        return 1;
    }
    
    std::ifstream in(argv[1], std::ifstream::in);
    std::vector<std::vector<double>> dataset;
    int cluster_cnt = std::atoi(argv[2]);
    //printf("%d\n", cluster_cnt);
    
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
    /*for (int i = 0; i < 5; i++) {
        for (double f : dataset[i]) {
            printf("%lf ", f);
        }
        printf("\n");
    }*/
 
    return 0;
}
