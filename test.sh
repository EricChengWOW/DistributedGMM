python3 createTestData.py
g++ -o seqGMM ./seqGMM.cpp
./seqGMM test_data.txt 5
python3 visualize.py
