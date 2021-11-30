#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <stdio.h>
#include <thrust/device_vector.h>
#include <thrust/system/system_error.h>
#include <vector>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <functional>
#include <map>
#include <omp.h>
#include <chrono>
#define CUDA_CALL(x) do { if((x) != cudaSuccess) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__); \
    printf("%s\n",cudaGetErrorString(x)); \
    system("pause"); \
    return EXIT_FAILURE;}} while(0)
static void HandleError(cudaError_t err, const char*file,int line){ 
if (err != cudaSuccess) {
 printf( "%s in %s at line %d\n", cudaGetErrorString( err ), file, line );
 exit( EXIT_FAILURE ); }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))
using namespace std;

/* funzione che sostituisce ogni valore del vettore di input con la loro posizione con ripetizione (stesso valore = stesso indice)
   per avere nodi con nomi da 1 ad n con n= numero di nodi del grafo */
vector<int> changeVect(vector<int> input, int N) {

	map<int, int> ranks; // creo un contenitore associativo chiave-valore per evitare di creare un vettore con dimensione = max valore nell'input

	int rank = 1;

	for (int index = 0; index < N; index++) { // riempio ranks associando a ogni valore la sua posizione

		int element = input[index];

		if (ranks[element] == 0) // aumento il rank solo se non si tratta di un valore ripetuto
		{
			ranks[element] = rank;
			rank++;
		}
	}

	for (int index = 0; index < N; index++) // sostituisco ogni valore del vettore di input con il relativo rank associato
	{
		int element = input[index];
		input[index] = ranks[input[index]];
	}

	return input;

}

/* funzione che apre il file di input .txt in formato SNAP e ricava i vettori IC e IR necessari per il calcolo dell'ExF e restituisce un intero che è il grado max nella rete*/
int calcoloICIR(string filename, vector <int>& IC, vector <int>& IR) {

	ifstream grafo_snap;  // apro il file contenente il grafo in formato SNAP
	grafo_snap.open(filename, ios::in);

	vector<int> c1;  // creo i vettori che conterranno i dati
	vector<int> c2;

	vector< pair <int, int> > vect; // creo vettore di coppie

	// controllo che il file esista
	if (!grafo_snap) {
		cout << "File non esistente";
	}
	else {

		int nodo_partenza; // dichiaro due int di supporto equivalenti a testa e coda dell'arco
		int nodo_arrivo;

		while (grafo_snap >> nodo_partenza && grafo_snap >> nodo_arrivo) {  // salvo il primo e il secondo elemento di ogni riga nei vettori
			vect.push_back(make_pair(nodo_partenza, nodo_arrivo));

			if (grafo_snap.eof()) // controllo che il ciclo scandisca il file fino al termine
				break;
		}
	}
	grafo_snap.close();

	sort(vect.begin(), vect.end());

	for (auto it = std::make_move_iterator(vect.begin()),
		end = std::make_move_iterator(vect.end()); it != end; ++it)
	{
		c1.push_back(std::move(it->first));
		c2.push_back(std::move(it->second));
	}

	vector<int> vect_union; // dichiaro il vettore unione per poter calcore il "rank"
	vect_union.reserve(c1.size() + c2.size()); // alloco la memoria per il vettore unione
	vect_union.insert(vect_union.end(), c1.begin(), c1.end()); // faccio l'unione dei due vettori
	vect_union.insert(vect_union.end(), c2.begin(), c2.end());

	vector<int> output; // vettore di output da dividire in seguito

	size_t N = vect_union.size();

	output = changeVect(vect_union, N);

	// devo calcolare IC e IR finali per darli in input alla funzione che calcola l'ExF

	// IC è già nella sua forma finale (ovvero = seconda meta' del vettore output), devo solo aggiungere uno 0 in testa
	for (int k = output.size() / 2; k < output.size(); k++) {
		IC.push_back(output[k]);
	}
	IC.insert(IC.begin(), 0);

	// I step per calcolo IR: calcolare il numero di nodi del grafo
	std::vector<int>::iterator it_numNodi;
	it_numNodi = std::max_element(output.begin(), output.end());
	int num_nodi = it_numNodi[0];

	// II step per calcolo IR: considero la prima parte del vettore output, da qui il nome preIR
	vector<int> preIR(output.begin(), output.begin() + output.size() / 2);

	// III step per calcolo IR: calcolo per ogni nodo il numero di ripetizioni (se ci sono)

	int j = 0;
	int i = 0;
	int sum = 1; // parte da 1 perche' ci sara' sempre almeno un nodo 
	int partialSum = 1;
	int maxSum = 1;

	while (j < preIR.size() && i < preIR.size() - 1) {


		if (preIR[i] != preIR[i + 1]) {

			IR.push_back(sum);
			j++;
			maxSum = max(partialSum, maxSum);
			partialSum = 1;
		}
		else {
			partialSum = partialSum + 1;

		}


		sum++;
		i++;



	}


	for (int k = j; k < num_nodi; k++) {
		IR.push_back(sum);
	}

	IR.insert(IR.begin(), 0); // inserisco in testa lo 0 come da definizione
	return  max(partialSum, maxSum);

}


/*
* funzione che calcola il grado dei cluster: ogni cluster è formato da nodo seme e altri due nodi (in questo caso non importa se i nodi sono a distanza 1 o a distanza 2)
* grado cluster=numero di archi che connettono nodi nel cluster con nodi non nel cluster
*/
__device__ int gradoCluster(int* nodi, int* IR, int* IC) {
	// per ogni nodo calcolo il numero di vicini
	int grado = 0;
	for (int a = 0; a < 3; a++) {
		int nodo = nodi[a];//nodo considerato per calcolare il numero di vicini
		int valRiga = IR[nodo] - IR[nodo - 1];//numero di vicini del nodo
		if (valRiga > 0) {//controllo se il numero di vicini è >0 altrimenti il grado è 0
			// printf("Numero di vicini %i:%i\n", nodo, valRiga);
			for (int k = IR[nodo - 1] + 1; k < IR[nodo - 1] + 1 + valRiga; k++) {//scandisco tutti i vicini del nodo considerato e controllo se questi nodi non sono gi? nel cluster
				int vicino = IC[k];
				// printf("Nodo vicino di %i:%i\n", nodo ,vicino);
				if (vicino != nodi[0] && vicino != nodi[1] && vicino != nodi[2]) {// se il vicino del nodo non è nel cluster aumento il grado del cluster
					grado++;
				}

			}
		}
	}

	//printf("Grado del cluster(%i, %i, %i):%i,\n", nodi[0], nodi[1], nodi[2], grado);

	return grado;
}

/*
* kernel che dati IC, IR calcola ExF per ogni nodo del grafo
*/
__global__ void expectedForce(int* IR_vec, int* IC_vec, int n_IR, int gradoMax, double* d_exf)
{

	//indice i per scorrere l'array IR_vec 
	int seed = blockDim.x * blockIdx.x + threadIdx.x;
	if (seed < n_IR && seed != 0) {//controllo di non eccedere la dimensione dell'array e che non venga calcolato ExF del primo elemento di IR_vec
		double ExF = 0; //ogni thread memorizza il valore di ExF in questa variabile
		int indiceGradi = 0;
		int totalFI = 0;//ogni thread memorizza la somma dei gradi in questa variabile
		//printf("Thread %d calcola nodo %d \n", seed, IR_vec[seed]);

		int valRiga = IR_vec[seed];
		int* distOne = new int[gradoMax];//vettore di vicini a distanza 1 dal nodo seed di dimensione gradoMax (numero massimo di vicini di un nodo nel grafo considerato)
		if (distOne == NULL) { printf("distOne failed\n"); return; } // controllo se il puntatore è nullo

		int dist = IR_vec[seed] - IR_vec[seed - 1]; //numero di vicini del nodo seed
		if (dist == 0) { //se il nodo considerato non ha vicini allora ExF=0
			d_exf[seed] = 0;
			//printf("Exf del nodo %i: %f\n", seed, d_exf[seed]);

		}
		else // se il nodo considerato ha almeno un vicino
		{

			//calcolo i nodi a distanza uno dal nodo seme e li metto nel vettore distOne
			int indiceDistOne = 0;
			for (int k = IR_vec[seed - 1] + 1; k <= valRiga; k++) {
				int valB = IC_vec[k];
				//printf("Nodo a distanza uno da %i:%i\n",seed, valB);
				if (valB != seed) {//controllo che il vicino del nodo non sia il nodo stesso 
					distOne[indiceDistOne] = valB;
					indiceDistOne++;
				}
			}
			if (indiceDistOne == 0) {//se a questo punto il vettore distOne è vuoto allora il nodo seme punta solo a se stesso


				d_exf[seed] = 0;
				//printf("Exf del nodo %i: %f\n", seed, d_exf[seed]);

			}
			else {//se il vettore distOne non è vuoto
				//calcolo la dimensione dell'array gradi considerando il caso peggiore: ogni nodo nodo ha gradoMax vicini e tutti i suoi vicini hanno gradoMax vicibi
				int worstCaseSize = gradoMax * gradoMax + gradoMax * (gradoMax - 1);
				int* gradi = new int[worstCaseSize]; //this allocates memory on a local memory runtime heap which has the lifetime of the context
				if (gradi == NULL) { printf("gradi failed\n"); return; }// controllo se il puntatore è nullo

				//devo creare i cluster con 2 nodi a distanza 1 da seed
			// il doppio ciclo for mi permette di trovare le combinazioni di nodi a distanza 1 da seed
				for (int k = 0; k < indiceDistOne - 1; k++) {
					for (int j = k + 1; j < indiceDistOne; j++) {
						int nodiCluster[3];  //vettore che contiene i 3 nodi che formano il cluster considerato
						//printf("Nodi nel cluster:%i,%i,%i\n", seed, distOne[k], distOne[j]);
						nodiCluster[0] = seed;
						nodiCluster[1] = distOne[k];
						nodiCluster[2] = distOne[j];

						int grado = gradoCluster(nodiCluster, IR_vec, IC_vec); //calcolo il grado del cluster considerato
						int mult = 2; //ogni cluster si puo' creare in due combinazioni 

						for (int count = 0; count < mult; count++) { // inserisco nel vettore gradi 2 volte il grado del cluster considerato
							gradi[indiceGradi] = grado;

							indiceGradi++;
						}


						totalFI += mult * grado;

					}
				}


				//devo creare i cluster con un nodo a distanza 1 e 1 a distanza 2 dal nodo seed
				for (int k = 0; k < indiceDistOne; k++) { //scandisco i nodi a distanza 1 dal nodo seed

					int nodo = distOne[k]; //considero un nodo a distanza 1 
					// printf("Secondo nodo cluster:%i\n", nodo);
					 //itero su tutti i nodi a distanza 1 e trovo i loro vicini
					int valRiga = IR_vec[nodo] - IR_vec[nodo - 1]; //numero di vicini del nodo considerato
					if (valRiga > 0) { //se il nodo considerato ha almeno un vicino posso procedere con il calcolo del cluster
						//printf("Numero di vicini di %i:%i\n", nodo, valRiga);
						 //printf("Inizio del for %i ", IR[nodo - 1] + 1);
						for (int j = IR_vec[nodo - 1] + 1; j < IR_vec[nodo - 1] + 1 + valRiga; j++) { //scandisco i vicini del nodo considerato che sono i nodi a distanza 2 del nodo seed
							//nodo vicino del vicino, calcolo il numero di archi uscenti 
							int valB = IC_vec[j];
							int nodiCluster[3];
							// printf("Terzo nodo cluster:%i\n", valB);

							// printf("Nodi nel cluster a dist 1 e 2:%i,%i,%i\n", seed, nodo, valB);
							nodiCluster[0] = seed;
							nodiCluster[1] = nodo;
							nodiCluster[2] = valB;;


							int grado = gradoCluster(nodiCluster, IR_vec, IC_vec); // calcolo il grado del cluster
							gradi[indiceGradi] = grado; //inserisco il grado calcolato nel vettore gradi


							indiceGradi++;

							totalFI += grado;
						}

					}


				}
				delete[]distOne;

				/*
				* Dopo aver calcolato il vettore dei gradi dei cluster controllo :
				* - se la somma totale dei gradi è 0 allora anche ExF=0
				* - altrimenti normalizzo il vettore del gradi (dividendo ogni grado per la somma dei gradi) e calcolo ExF
				*/

				if (totalFI == 0) {

					d_exf[seed] = 0;
					//printf("Exf del nodo %i: %f\n", seed, d_exf[seed]);

				}
				else {

					double norm = 0;

					for (int K = 0; K < indiceGradi; K++) { //normalizzo l'array gradi e calcolo il valore di ExF
						if (gradi[K] != 0) {
							norm = (float)gradi[K] / totalFI;
							ExF -= log(norm) * norm;
						}


					}


					d_exf[seed] = ExF;
					//printf("Exf del nodo %i: %f\n", seed, d_exf[seed]);

				}

				delete[] gradi;

			}

		}


	}



}


int main(int argc, char* argv[])
{
	auto start = chrono::steady_clock::now();
	//definisco le variabili per calcolare l'occupazione della memoria
	size_t free_byte;
	size_t total_byte;

	//la variabile malloc_limit viene inizializzata alla quantità di memoria heap necessaria per il calcolo di ExF 
	//(la memoria heap ci serve per memorizzare il vettore dei vicini a distanza 1 e il vettore dei gradi dei cluster (uno per ogni thread))
	size_t malloc_limit = 1073741824; // 512 MB = 536870912; 1 GB = 1073741824; 256 MB = 268435456;
	HANDLE_ERROR(cudaDeviceSetLimit(cudaLimitMallocHeapSize, malloc_limit)); //imposto la capienza della memoria heap a mallo_limit
	//printf("dimensione della heap:%zi\n", malloc_limit);

	//creo i vettori h_IC e h_IR nell'host
	vector <int> h_IC;
	vector <int> h_IR;
	string name_file(argv[1]);
	//inizializzo i vettori h_IC e h_IR in base al file .txt di input tramite la funzione calcoloICIR 
	//la funzione calcoloICIR restituisce maxDegree che è un valore che verrà riutilizzato nel kernel per trovare un upper bound alla dimensione degli array distOne
	int maxDegree = calcoloICIR(name_file,h_IC, h_IR);
	//printf("Max degree:%i\n", maxDegree);

	/*cout << "IR = ";
	cout << "[";
	for (int i = 0; i < h_IR.size() - 1; i++) {
		cout << h_IR[i] << ", ";
	}
	cout << h_IR[h_IR.size() - 1] << "]" << endl;


	cout << "IC = ";
	cout << "[";
	for (int i = 0; i < h_IC.size() - 1; i++)
	{
		cout << h_IC[i] << ", ";
	}
	cout << h_IC[h_IC.size() - 1] << "]";*/

	//calcolo la lunghezza dei vettori h_IC e h_IR
	size_t n_IC = h_IC.size();
	size_t n_IR = h_IR.size();


	//per poter passare i vettori h_IC e h_IR al kernel devo creare delle copie (d_IC e d_IR) del tipo device_vector della libreria thrust
	thrust::device_vector<int> d_IR = h_IR;
	thrust::device_vector<int> d_IC = h_IC;

	//creo i puntatori raw che verranno passati al kernel per i vettori h_IC e h_IR 
	int* IR_vec = thrust::raw_pointer_cast(d_IR.data());
	int* IC_vec = thrust::raw_pointer_cast(d_IC.data());


	int N = n_IR;
	//printf("Numero di nodi del grafo:%i\n", N-1);

	size_t size = N * sizeof(double);

	//Creo il vettore h_exf in cui verranno memorizzati i valori ExF di tutti i nodi del grafo e lo alloco nella memoria dell'host
	double* h_exf = (double*)malloc(size);
	memset(h_exf, 0, N);//imposto a 0 gli elementi di h_exf

	//Creo il vettore d_exf che verrà usato dal device per memorizzare i valori di ExF di tutti i nodi del grafo, lo alloco in memoria host e lo copio in memoria device
	double* d_exf;
	HANDLE_ERROR(cudaMalloc(&d_exf, size));
	HANDLE_ERROR(cudaMemcpy(d_exf, h_exf, size, cudaMemcpyHostToDevice));

	//lancio il kernel expectedForce
	int threadsPerBlock = 512;
	int blocksPerGrid = (n_IR + threadsPerBlock - 1) / threadsPerBlock;
	expectedForce << <blocksPerGrid, threadsPerBlock >> > (IR_vec, IC_vec, n_IR, maxDegree, d_exf);
	CUDA_CALL(cudaDeviceSynchronize());

	//copio i risultati da d_exf a h_exf e li stampo 
	HANDLE_ERROR(cudaMemcpy(h_exf, d_exf, size, cudaMemcpyDeviceToHost));
	for (int i = 1; i < n_IR; i++) {
		printf("Exf del nodo %i: %f\n", i, h_exf[i]);
	}

	//libero la memoria allocata
	cudaFree(d_exf);
	free(h_exf);

	// show memory usage of GPU
	HANDLE_ERROR(cudaMemGetInfo(&free_byte, &total_byte));

	double free_db = (double)free_byte;

	double total_db = (double)total_byte;

	double used_db = total_db - free_db;

	printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n",

		used_db / 1024.0 / 1024.0, free_db / 1024.0 / 1024.0, total_db / 1024.0 / 1024.0);



	auto end = chrono::steady_clock::now();
	auto diff = end - start;//calcolo tempo totale di computazione
	cout << chrono::duration <double, milli>(diff).count() << " ms" << endl;

	return 0;
}
