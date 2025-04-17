#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
using namespace std;
int main() {
	string filename;
    vector<double> x, y, sigmay;
	double xval,yval,sigma;
	int scelta;
	
	cout<<"___________REGRESSIONI LINEARI PESATE______ "<<endl;
    cout<<"Seleziona se i dati possiedono una particolare distribuzione di probabilita'': "<<endl;

    cout<<"Opzioni: - inserisci 0 per inserire le incertezze delle y personalizzate        - Inserisci 1 per dati con distribuzione di Possion "<<endl;
  	cin>>scelta;
		
    //operazioni dal file e presa dati
	cout << "Inserire il nome del file di input: ";
    cin >> filename;
    ifstream fin(filename);
    if(!fin){
    	cout<<"Il file non è stato caricato correttamente"<<endl;
    	return -1;
	}
	
	if(scelta==0){
		while(fin>>xval>>yval>>sigma){
		if(xval>=320&&xval<=335 || xval>=390&&xval<=405){
			x.push_back(xval);
    			y.push_back(yval);
    			sigmay.push_back(sigma);	
		}
		}	
	}
	if(scelta==1){
		while(fin>>xval>>yval){
		if(xval>=320&&xval<=335 || xval>=390&&xval<=405){
			x.push_back(xval);
    			y.push_back(yval);
		}		
		}
		for(int i=0;i<y.size();i++){
			sigmay.push_back(sqrt(y.at(i)));
		}
	}
   
	    
    //Analisi dati         
	// Calcolo delle sommatorie utili
    double spesi = 0.0;   
    double sx = 0.0;   
    double sy = 0.0;  
    double sxx = 0.0;  
    double sxy = 0.0;  
	
	
	//calcolo somme intermedie
    for (int i = 0; i < x.size(); i++) {
        double w = 1.0 / (sigmay[i] * sigmay[i]); // peso = 1/s_i^2
        spesi  += w; //sommatoria pesi 
        sx  += x.at(i) * w; 
        sy  += y.at(i)* w;
        sxx += x.at(i) * x.at(i) * w;
        sxy += x.at(i) * y.at(i) * w;
    }

    // Calcolo del DELTA
    double delta = (spesi * sxx) - (sx * sx);

    // Stima dei parametri a e b
    double a = ( (sxx * sy) - (sx * sxy) ) / delta;
    double b = ( (spesi  * sxy) - (sx * sy ) ) / delta;

    //Calcolo delle incertezze relative ai parametri   
    double sigma_a = sqrt(sxx/ delta);
    double sigma_b = sqrt(spesi/ delta);

    // Stampa dei risultati
    cout << "=== RISULTATI DEL FIT LINEARE PESATO per una retta y=a+bx ===" << endl;
    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "sigma_a = " << sigma_a <<endl;
    cout << "sigma_b = " << sigma_b << endl;
	
    string scelta2;
    cout<<"I due parametri a e b sono covarianti?: ( scrivi: si oppure no )  "<<endl;
    cin>>scelta2;

    return 0;
}

