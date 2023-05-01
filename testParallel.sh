python3 createTestData.py
make clean && make
mpirun -n 2 ./parallelGMM test_data.txt 5
python3 visualize.py
