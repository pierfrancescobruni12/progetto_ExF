
#!/bin/sh
echo "Codice seriale"
g++ -std=c++11 ./calcoloExF_seriale.cpp -o seriale
./seriale grafiSNAP/grafo1.txt

echo "Codice parallelo_1"
nvcc 1_calcoloExF_parallelo.cu -o exf_parallelo_1
./exf_parallelo_1 grafiSNAP/grafo1.txt

echo "Codice parallelo_2"
nvcc -rdc=true 2_calcoloExF_parallelo.cu -lcudadevrt -o exf_parallelo_2
./exf_parallelo_2 grafiSNAP/grafo1.txt

echo "Codice parallelo_3"
nvcc -rdc=true 3_calcoloExF_parallelo.cu -lcudadevrt -o exf_parallelo_3
./exf_parallelo_3 grafiSNAP/grafo1.txt