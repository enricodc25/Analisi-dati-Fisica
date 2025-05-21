#include<iostream>
#include<vector>
#include<iomanip>
#include<cmath>
#include<fstream>
#include <cstdlib>
#include<string>
using namespace std;

double deviazioneStandardCampionaria(vector<double> dati);
double media(vector<double> v);

void errore_nonstatiticità(vector<double> tdil,vector<double> tcomp,double bdil, double bcomp){
    
    double tcompmedia=media(tcomp);
    double tdilmedia=media(tdil);
    double r,incertezza_n;
    double median;
    double ndil,ncomp;
    double tmedia,errorereltemp;

    cout<<"DEVO CALCOLARE L'ERRORE DELLA NON STATICITA, HO BISOGNO DEL VALORE DI R "<<endl;
    r=8.314;
    //Ora la temperatura la converto in kelvin
    ndil=(bdil/(r*(tdilmedia+273.15)))*9.806/100;
    ncomp=(bcomp/(r*(tcompmedia+273.15)))*9.806/100;
    //per la n media caratteristica di una temperatura faccio la media tra le moli di compressione e di dilatazione
    median=(ndil+ncomp)/2;

    incertezza_n=abs((ndil-ncomp)/((2*sqrt(3))*median));
    //per la temperatura media faccio la media tra le due medie della temperatura di compressione e di dilatazione
    tmedia=(tdilmedia+tcompmedia)/2;

    for(int i=0;i<tdil.size();i++){
        tcomp.push_back(tdil.at(i));
    }
    double devstandardtemp=deviazioneStandardCampionaria(tcomp);
    double erroredevt=0.3/(sqrt(3)*tmedia);
    
    cout<<"La temperatura media per questa particolare temperatura: "<<tmedia<<endl;
    cout<<"Variazione relativa della temperatura:  "<<devstandardtemp/tmedia<<endl;
    cout<<"Deviazione standard della temperatura sapendo l'errore massimo sulla temperatura: "<<erroredevt<<endl;
    cout<<"Moli in compressione: "<<ncomp<<endl;
    cout<<"Moli in dilatazione: "<<ndil<<endl;
    cout<<"Numero di moli medio tra compressione e dilatazione: "<<median<<endl;
    cout<<"Incertezza relativa su n dovuta alla differenza tra dilatazione e compressione: "<<incertezza_n<<endl;
}

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

double media(vector<double> v){
    double somma=0.0;
    for(int i=0;i<v.size();i++){
        somma+=somma+v.at(i);
    }

    return somma/v.size();
}


int main(){
    string f1,f2;
    double bcomp,bdil;
    vector<double> x1,y1,z1;
    vector<double> x2,y2,z2;
    double xval,yval,zval;


    cout<<"*****CALCOLO INFORMAZIONI MOLE AD UNA TEMPERATURA PRECISA***"<<endl;
    cout<<"Inserisci il nome del file della compressione: ";
    cin>>f1;
    cout<<"inserisci il nome del file della dilatazione: ";
    cin>>f2;

    ifstream fin1(f1);
    ifstream fin2(f2);
    cout<<"Inserisci la b della compressione: ";
    cin>>bcomp;
    cout<<"Inserisci la b della dilatazione: ";
    cin>>bdil;

    while(fin1>>xval>>yval>>zval){
	    x1.push_back(xval);    		
    	y1.push_back(yval);
    	z1.push_back(zval);
    
    }  

    while(fin2>>xval>>yval>>zval){
	    x2.push_back(xval);    		
    	y2.push_back(yval);
    	z2.push_back(zval);
    
    }     
    errore_nonstatiticità(z2,z1,bdil,bcomp);

    return 0;
}
