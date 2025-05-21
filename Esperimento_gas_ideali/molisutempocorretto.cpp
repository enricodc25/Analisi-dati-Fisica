#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
using namespace std;
void interpolazionesemplice(vector<double> x, vector<double> y, double &a, double &b);
double deviazioneStandardCampionaria(vector<double> dati);

int main() {

    std::ifstream input("prova_tenuta.txt");
    std::ofstream output("prova.txt");
    vector<double> moli,tempi;
    double n;
    double alfa,nzero;
    if (!input) {
        std::cerr << "Errore nell'apertura del file dati.txt\n";
        return 1;
    }

    if (!output) {
        std::cerr << "Errore nella creazione del file risultato.txt\n";
        return 1;
    }

    // Costanti note
    const double R = 8.314;      // [J/(mol·K)]
    const double V_0 = 0.4475;    // [m^3], esempio

    double p, V_CL, T;
    int i = 0;
    double t_i;

    while (input >> p >> V_CL >> T) {
        t_i = 0.1 * i;
        n = (p* (V_CL + V_0)*9.806) / (R * (T+273.15)*100);
        moli.push_back(n);
        tempi.push_back(t_i);
        output << t_i << "\t" << n << "\n";
        i++;
    }

    input.close();
    output.close();

    std::cout << "Calcolo completato. File risultato.txt generato.\n";
    interpolazionesemplice(tempi,moli,nzero,alfa);
    cout<<"---CALCOLO ALFA, errore sistematico dovuto ad una perdita di gas lineare nel tempo"<<endl;
    cout<<" Alfa: "<<alfa<<"  Nzero: "<<nzero<<endl;
    return 0;
}
//INTERPOLAZIONE SEMPLICE CON INCERTEZZE SULLE Y TUTTE UGUALI
void interpolazionesemplice(vector<double> x, vector<double> y, double &a, double &b){  
    double sx = 0.0;   
    double sy = 0.0;  
    double sxx = 0.0;  
    double sxy = 0.0; 
    double n=x.size(); 
    double sigmay=deviazioneStandardCampionaria(y);

	//calcolo somme intermedie
    for (int i = 0; i < n; i++) {
        sx  += x.at(i); 
        sy  += y.at(i);
        sxx += x.at(i) * x.at(i);
        sxy += x.at(i) * y.at(i);
    }

    // Calcolo del DELTA
    double delta = n*(sxx) - (sx * sx);

    // Stima dei parametri a e b
    a = ( (sxx * sy) - (sx * sxy) ) / delta;
    b = ( (n* sxy) - (sx * sy ) ) / delta;

    //Calcolo delle incertezze relative ai parametri   
    double sigma_a = sigmay*sqrt(sxx/ delta);
    double sigma_b = sigmay*sqrt(n/ delta);
}

// Deviazione standard campionaria
double deviazioneStandardCampionaria(vector<double> dati) {
    int n = dati.size();
    double somma = 0.0;
    for (double valore : dati) {
        somma += valore;
    }
    double media = somma / n;
    double sommaQuadrati = 0.0;

    for (double valore : dati) {
        sommaQuadrati += (valore - (somma/n)) * (valore - (somma/n));
    }
    return sqrt(sommaQuadrati / (n - 1));
}