#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>

using namespace std;

int main() {
    vector<double> A_values, S_values, p, q;

    // --- 1. LETTURA DATI ---
    ifstream inputFile("datipqcorretti.txt");
    if (!inputFile) {
        cerr << "Errore: Impossibile aprire il file datipqcorretti.txt" << endl;
        return 1;
    }

    double valueP, valueQ;
    while (inputFile >> valueP >> valueQ) {
        p.push_back(valueP);
        q.push_back(valueQ);
    }
    inputFile.close();

    // Controllo lettura dati
    cout << "Dati letti correttamente. Numero punti: " << p.size() << endl;

    // --- 2. PREPARAZIONE SCANSIONE (GRID SEARCH) ---
    ofstream outputFile("residui_A.dat");
    double A_center = 147.941; //puoi divertirti a modificarlo per vedere come cambia il minimo, indica i valori visuallizzati sull'asse x del grafico

    double A_min_scan = A_center - 2;
    double A_max_scan = A_center + 2;
    double step = 0.001;
    
    double S_min = numeric_limits<double>::max();
    double A_best = 0.0;

    // --- 3. CICLO DI MINIMIZZAZIONE ---
    // Scansioniamo i valori di A per trovare quello che rende minima la somma S
    for (double A = A_min_scan; A <= A_max_scan; A += step) {
        double S = 0.0;

        for (int i = 0; i < p.size(); i++) {
            // Evita divisioni per zero nel modello q = (A*p)/(p-A)
            if (fabs(p[i] - A) < 1e-6) continue;

            // Calcolo del residuo: differenza tra dato sperimentale e modello
            double residuo = q[i] - (A * p[i]) / (p[i] - A);
            S += residuo * residuo; 
        }

        // Memorizziamo i risultati per l'analisi successiva
        A_values.push_back(A);
        S_values.push_back(S);
        outputFile << A << " " << S << endl;

        // Identificazione del minimo assoluto
        if (S < S_min) {
            S_min = S;
            A_best = A;
        }
    }
    outputFile.close();

    // --- 4. RICERCA INTERVALLO DI CONFIDENZA ---
    // Cerchiamo i valori di A dove S(A) = S_min + 1 (corrisponde a 1 sigma)
    double Ameno1 = 0, Apiu1 = 0;
    bool trovato_inf = false;

    for (size_t i = 0; i < A_values.size(); ++i) {
        // Nota: Il confronto diretto == con double è rischioso per la discretizzazione.
        // Idealmente si cerca il valore più vicino o si usa un'interpolazione.
        if (!trovato_inf && S_values[i] <= S_min + 1.0) {
            Ameno1 = A_values[i];
            trovato_inf = true;
        }
        if (trovato_inf && S_values[i] > S_min + 1.0 && A_values[i] > A_best) {
            Apiu1 = A_values[i];
            break;
        }
    }

    // Sovrascrittura manuale (come da tuo codice originale per precisione)
Ameno1 = 147.848; //questi valori li ho messi io a mano, per velocità
    Apiu1 = 148.034;
    double incertezza = (Apiu1 - Ameno1) / 2.0;

    // --- 5. OUTPUT RISULTATI ---
    cout << "\nValore di A che minimizza la somma dei residui:" << endl;
    cout << "A_best = " << A_best << endl;
    cout << "S_min  = " << S_min << endl;
    cout << "Incertezza su A_best: +/- " << incertezza << endl;
    cout << "Intervallo di confidenza (S_min + 1): [" << Ameno1 << ", " << Apiu1 << "]" << endl;

    return 0;
}