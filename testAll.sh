python3 createTestData.py
make clean && make
./seqGMM test_data.txt 5
mpirun -n 1 ./parallelGMM test_data.txt 5
mpirun -n 2 ./parallelGMM test_data.txt 5
mpirun -n 4 ./parallelGMM test_data.txt 5
mpirun -n 8 ./parallelGMM test_data.txt 5
mpirun -n 16 ./parallelGMM test_data.txt 5
mpirun -n 32 ./parallelGMM test_data.txt 5
mpirun -n 64 ./parallelGMM test_data.txt 5
mpirun -n 128 ./parallelGMM test_data.txt 5
# python3 visualize.py
