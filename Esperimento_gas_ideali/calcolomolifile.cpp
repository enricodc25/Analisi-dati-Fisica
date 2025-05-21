#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
using namespace std;
void interpolazionesemplice(vector<double> x, vector<double> y, double &a, double &b);
double deviazioneStandardCampionaria(vector<double> dati);
double media(vector <double> x);


int main() {

    vector<double> molicomp,molidil,tempicomp,tempidil;
    double n;
    string nome1,nome2;
    double v1,v2;
    const double R = 8.314;      // [J/(mol·K)]
    double p, V_CL, T;
    int i = 0;
    double t_i,nzero1,nzero2,alfa1,alfa2;
    double tvecchio;
    double alfa;
    
//ATTENZIONE CHE L'ALFA SI FA SUL FILE DI TENUTA




    cout<<"Inserisci il nome del file della compressione: ";
    cin>>nome1;
    cout<<"Inserisci il nome del file della dilatazione: ";
    cin>>nome2;

    ifstream fincomp(nome1);
    ifstream findilat(nome2);

    cout<<"Inserisci il volume morto compressione: ";
    cin>>v1;
    cout<<"Inserisci il volume morto dilatazione: ";
    cin>>v2;

    cout<<"Inserisci il tempo medio vecchio ( se la temperatura analizzata è quella di 0 gradi, inserire 0) ";
    cin>>tvecchio;
    // cout<<"Inserisci il valore di alfa: ";
    // cin>>alfa;
    alfa=-0.00000000336903;

    //carico ed elaboro i dati della compressione
    while (fincomp>> p >> V_CL >> T) {
        t_i = 0.1 * i;
        n = (p* (V_CL + v1)*9.806) / (R * (T+273.15)*100);
        molicomp.push_back(n);
        tempicomp.push_back(t_i);
        i++;
    }
    double j=0;
    //carico ed elaboro i dati della dilatazione
    while (findilat >> p >> V_CL >> T) {
        t_i = 0.1 * j;
        n = (p* (V_CL + v2)*9.806) / (R * (T+273.15)*100);
        molidil.push_back(n);
        tempidil.push_back(t_i);
        j++;
    }
    double molic=media(molicomp);
    double molid=media(molidil);
    
    //metto insieme i tempi della compressione e dilatazione
    // for(int i=0; i<tempidil.size();i++){
    //     tempicomp.push_back(tempidil.at(i));
    // }

    // //faccio la media delle moli nel caso di una temperatura precisa:
    // for(int i=0; i<molidil.size();i++){
    //     molicomp.push_back(molidil.at(i));
    // }
    double mediamoli=media(molicomp);

    double errorerelativomoli_comp,errorerelativomoli_dil;
    errorerelativomoli_comp=-(alfa/molic)*(tvecchio+(tempicomp.at(tempicomp.size()-1))/2);
    errorerelativomoli_dil=-(alfa/molid)*(tvecchio+tempicomp.at(tempicomp.size()-1)+((+tempidil.at(tempidil.size()-1)))/2);

    cout<<endl;
    cout<<"***************************************RISULTATI******************************"<<endl<<endl<<endl;
    cout<<"---Tempo totale compressione + dilatazione (diventera' il tempo vecchio per la temperatura successiva): " 
    <<(tempicomp.at(tempicomp.size()-1))+(tempidil.at(tempidil.size()-1))<<endl<<endl;
    cout<<"--- Tempo totale dall'inizio della presa dati: "<<tvecchio+(tempicomp.at(tempicomp.size()-1))+(tempidil.at(tempidil.size()-1))<<endl<<endl<<endl;
    cout<<"_______COMPRESSIONE_____"<<endl;
    cout<<"Moli medie: "<<molic<<endl;
    cout<<"Errore relativo delle moli : "<<errorerelativomoli_comp<<endl;
    cout<<"_______DILATAZIONE_____"<<endl;
    cout<<"Moli medie: "<<molid<<endl;
    cout<<"Errore relativo delle moli : "<<errorerelativomoli_dil<<endl<<endl;

    return 0;
}


double media(vector<double> x){
    double media;
    double somma=0;
    for(int i=0; i<x.size(); i++){
        somma+=x.at(i);
    }
    media=somma/x.size();

    return media;
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