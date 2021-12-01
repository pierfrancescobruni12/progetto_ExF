#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <functional>
#include <map>
#include <chrono>
#include <cstdio>




using namespace std;

/* funzione che sostituisce ogni valore del vettore di input con la sua posizione con ripetizione (stesso valore = stesso indice)
   per avere nodi con nomi da 1 ad n con n = numero di nodi del grafo */
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

/* funzione che apre il file di input .txt in formato SNAP e ricava i vettori IC e IR necessari per il calcolo dell'ExF */
void calcoloICIR(string filename, vector <int>& IC, vector <int>& IR) {

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

	int N = vect_union.size();

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

	while (j < preIR.size() && i < preIR.size() - 1) {

		if (preIR[i] != preIR[i + 1]) {
			IR.push_back(sum);
			j++;
		}
		sum++;
		i++;
	}

	for (int k = j; k < num_nodi; k++) {
		IR.push_back(sum);
	}

	IR.insert(IR.begin(), 0); // inserisco in testa lo 0 come da definizione
}


/*
* funzione che calcola il grado dei cluster: ogni cluster è formato da nodo seme e altri due nodi (in questo caso non importa se i nodi sono a distanza 1 o a distanza 2)
* grado cluster=numero di archi che connettono nodi nel cluster con nodi non nel cluster
*/
int gradoCluster(vector<int> nodi, vector<int> IR, vector<int> IC) {
	// per ogni nodo calcolo il numero di vicini
	int grado = 0;
	for (int a = 0; a < nodi.size(); a++) {
		int nodo = nodi[a]; // nodo considerato per calcolare il numero di vicini
		int valRiga = IR[nodo] - IR[nodo - 1]; // numero di vicini del nodo
		if (valRiga > 0) { // controllo se il numero di vicini è >0 altrimenti il grado è 0
			// printf("Numero di vicini %i:%i\n", nodo, valRiga);
			for (int k = IR[nodo - 1] + 1; k < IR[nodo - 1] + 1 + valRiga; k++) { // scandisco tutti i vicini del nodo considerato e controllo se questi nodi non sono già nel cluster
				int vicino = IC[k];
				// printf("Nodo vicino di %i:%i\n", nodo ,vicino);
				if (vicino != nodi[0] && vicino != nodi[1] && vicino != nodi[2]) { // se il vicino del nodo non è nel cluster aumento il grado del cluster
					grado++;
				}

			}
		}
	}
	//printf("Grado del cluster con seed %i (%i, %i, %i):%i,\n", nodi[0],nodi[0], nodi[1], nodi[2], grado);

	return grado;
}

/*
* funzione che dati IC, IR e il nodo seed calcola ExF per il nodo seed
* questa funzione poi sarà chiamata per ogni nodo del grafo
*/
double expectedForce(vector<int> IR, vector<int> IC, int i) {

	vector<int> gradi; // vettore che conterrà i gradi del cluster considerando il nodo seed (i)
	double ExF(0);
	int totalFI = 0;
	int seed = i;
	int valRiga = IR[seed]; // considero valRiga come l'indice massimo da leggere in IC

	vector<int> distOne; // vettore che contiene i nodi a distanza 1 dal nodo seed
	int dist = IR[seed] - IR[seed - 1]; // numero di vicini del nodo seed


	if (dist == 0) { // se il nodo considerato non ha vicini allora ExF=0
		ExF = 0;
	}
	else // se il nodo considerato ha almeno un vicino
	{
		// calcolo i nodi a distanza uno dal nodo seme e li metto nel vettore distOne
		for (int k = IR[seed - 1] + 1; k <= valRiga; k++) {
			int valB = IC[k];
			//printf("Nodo a distanza uno da %i:%i\n", seed, valB);
			if (valB != seed) {//devo controllare i cicli
				distOne.push_back(valB);
			}
		}
		if (distOne.empty()) {//se il vettore distOne è vuoto (nel caso in cui un nodo punti solo a se stesso=)
			ExF = 0;
		}
		else {

			// devo creare i cluster con 2 nodi a distanza 1 da seed
			// il doppio ciclo for mi permette di trovare le combinazioni di nodi a distanza 1 da seed
			for (int k = 0; k < distOne.size() - 1; k++) {
				for (int j = k + 1; j < distOne.size(); j++) {
					vector<int> nodiCluster; // vettore che contiene i 3 nodi che formano il cluster considerato
					//printf("Nodi nel cluster:%i,%i,%i\n", seed, distOne[k], distOne[j]);
					nodiCluster.push_back(seed);
					nodiCluster.push_back(distOne[k]);
					nodiCluster.push_back(distOne[j]);

					int grado = gradoCluster(nodiCluster, IR, IC); // calcolo il grado del cluster considerato
					int mult = 2; // ogni cluster si puo' creare in due combinazioni 
					/*if (areConnected(nodiCluster, IR, IC)) {
						mult = 4;
					} */
					//printf("Combinazioni del cluster:%i\n", mult);
					for (int count = 0; count < mult; count++) { // inserisco nel vettore gradi 2 volte il grado del cluster considerato
						gradi.push_back(grado);

					}

					totalFI += mult * grado;

				}
			}

			// devo creare i cluster con un nodo a distanza 1 e 1 a distanza 2 dal nodo seed
			for (int k = 0; k < distOne.size(); k++) { //scandisco i nodi a distanza 1 dal nodo seed

				int nodo = distOne[k]; // considero un nodo a distanza 1 
				// printf("Secondo nodo cluster:%i\n", nodo);
				 // itero su tutti i nodi a distanza 1 e trovo i loro vicini
				int valRiga = IR[nodo] - IR[nodo - 1]; //numero di vicini del nodo considerato
				if (valRiga > 0) { //se il nodo considerato ha almeno un vicino posso procedere con il calcolo del cluster
					//printf("Numero di vicini di %i:%i\n", nodo, valRiga);
					 //printf("Inizio del for %i ", IR[nodo - 1] + 1);
					for (int j = IR[nodo - 1] + 1; j < IR[nodo - 1] + 1 + valRiga; j++) { // scandisco i vicini del nodo considerato che sono i nodi a distanza 2 del nodo seed
						// nodo vicino del vicino, calcolo il numero di archi uscenti 
						int valB = IC[j];
						vector<int> nodiCluster;
						// printf("Terzo nodo cluster:%i\n", valB);

						//printf("Nodi nel cluster a dist 1 e 2:%i,%i,%i\n", seed, nodo, valB);
						nodiCluster.push_back(seed);
						nodiCluster.push_back(nodo);
						nodiCluster.push_back(valB);

						int grado = gradoCluster(nodiCluster, IR, IC); // calcolo il grado del cluster
						gradi.push_back(grado); // inserisco il grado calcolato nel vettore gradi


						totalFI += grado;
					}

				}


			}

			/*
			* Dopo aver calcolato il vettore dei gradi dei cluster controllo :
			* - se la somma totale dei gradi è 0 allora anche ExF=0
			* - altrimenti normalizzo il vettore del gradi (dividendo ogni grado per la somma dei gradi) e calcolo ExF
			*/

			if (totalFI == 0) {
				ExF = 0;
			}
			else {
				double norm = 0;
				//printf("valori gradi per seed %i:", seed);
				for (int i = 0; i < gradi.size(); i++) {
					if (gradi[i] != 0) {

						norm = (float)gradi[i] / totalFI;
						ExF -= log(norm) * norm;
					}

				}
			}


		}
	}
	//printf("Numero di vicini:%i\n", distOne.size());
	//printf("Lunghezza vettore gradi:%i\n", gradi.size());

	return ExF;
}

int main(int argc, char* argv[]) {

	auto start = chrono::steady_clock::now();
	vector <int> IC;
	vector <int> IR;

	string name_file(argv[1]);
	//inizializzo i vettori h_IC e h_IR in base al file .txt di input tramite la funzione calcoloICIR 
	//la funzione calcoloICIR restituisce maxDegree che è un valore che verrà riutilizzato nel kernel per trovare un upper bound alla dimensione degli array distOne
	calcoloICIR(name_file,IC, IR);
	/*cout << "IR = ";
	cout << "[";
	for (int i = 0; i < IR.size() - 1; i++) {
		cout << IR[i] << ", ";
	}
	cout << IR[IR.size() - 1] << "]" << endl;


	cout << "IC = ";
	cout << "[";
	for (int i = 0; i < IR.size() - 1; i++)
	{
		cout << IC[i] << ", ";
	}
	cout << IC[IC.size() - 1] << "]";*/


	// Per ogni nodo del grafo letto calcolo ExF 

	double exf;
	for (int seed = 1; seed <= IR.size() - 1; seed++) {
		exf = expectedForce(IR, IC, seed);
		printf("Exf del nodo %i: %f\n", seed, exf);

	}

	auto end = chrono::steady_clock::now();
	auto diff = end - start;
	// cout << chrono::duration <double, milli>(diff).count() << " ms" << endl;


	return 0;
}
